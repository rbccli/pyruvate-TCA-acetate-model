library(cmdstanr)
library(posterior)

# data_list must already exist in the R environment.
# It should be the same data list used for the posterior Stan fit.

n_prior <- 5000L
n_chains <- 4L

stopifnot(n_prior %% n_chains == 0L)

prior_rng_model <- cmdstan_model(
  "tca_ferm_prior_rng_fixed_param.stan"
)

prediction_model <- cmdstan_model(
  "tca_ferm_prior_posterior_predictive.stan"
)

# Stage 1: draw independent samples directly from the priors.
# This uses no warmup and no NUTS, so it should finish almost immediately.
prior_parameter_draws <- prior_rng_model$sample(
  seed = 20260816,
  chains = n_chains,
  parallel_chains = n_chains,
  iter_warmup = 0,
  iter_sampling = n_prior %/% n_chains,
  fixed_param = TRUE,
  refresh = 50
)

# Stage 2: for each prior draw, run the original Stan ODE exactly once and
# generate y_rep and f_rep with the same observation model as the posterior.
prior_data <- data_list
prior_data$prior_only <- 1L

prior_gq <- prediction_model$generate_quantities(
  data = prior_data,
  fitted_params = prior_parameter_draws$output_files(),
  seed = 20260817,
  parallel_chains = n_chains
)

prior_y_rep <- posterior::as_draws_matrix(
  prior_gq$draws("y_rep")
)

prior_f_rep_all <- posterior::as_draws_matrix(
  prior_gq$draws("f_rep")
)

saveRDS(
  list(
    prior_y_rep = prior_y_rep,
    prior_f_rep_all = prior_f_rep_all
  ),
  file = "stan_prior_predictive_draws_fast.rds"
)

cat(
  "Generated",
  nrow(prior_y_rep),
  "prior predictive draws.\n"
)



library(posterior)

set.seed(20260921)

# ------------------------------------------------------------
# Extract existing posterior model predictions
# ------------------------------------------------------------

post_y_hat <- posterior::as_draws_matrix(
  tca_ferm_fit_real$draws("y_hat")
)

post_f_hat <- posterior::as_draws_matrix(
  tca_ferm_fit_real$draws("f_hat")
)

post_sigma_c <- as.numeric(
  posterior::as_draws_matrix(
    tca_ferm_fit_real$draws("sigma_c")
  )[, 1]
)

post_sigma_f <- as.numeric(
  posterior::as_draws_matrix(
    tca_ferm_fit_real$draws("sigma_f")
  )[, 1]
)


# Put y_hat and f_hat columns in numerical order

y_hat_names <- paste0(
  "y_hat[",
  data_list$N_t,
  ",",
  seq_len(data_list$NY),
  "]"
)

f_hat_names <- paste0(
  "f_hat[",
  1:14,
  "]"
)

stopifnot(
  all(y_hat_names %in% colnames(post_y_hat)),
  all(f_hat_names %in% colnames(post_f_hat))
)

post_y_hat <- post_y_hat[
  ,
  y_hat_names,
  drop = FALSE
]

post_f_hat <- post_f_hat[
  ,
  f_hat_names,
  drop = FALSE
]


# Use 5000 posterior draws for faster density plotting

n_post_plot <- min(
  5000,
  nrow(post_y_hat)
)

post_keep <- sample.int(
  nrow(post_y_hat),
  n_post_plot
)

post_y_hat <- post_y_hat[post_keep, , drop = FALSE]
post_f_hat <- post_f_hat[post_keep, , drop = FALSE]
post_sigma_c <- post_sigma_c[post_keep]
post_sigma_f <- post_sigma_f[post_keep]


# ------------------------------------------------------------
# Generate posterior predictive replicates
# using the same log10 observation models as the likelihood
# ------------------------------------------------------------

posterior_y_rep <- vapply(
  seq_len(data_list$NY),
  function(j) {
    10^rnorm(
      n_post_plot,
      mean = log10(
        pmax(post_y_hat[, j], 1e-12)
      ),
      sd = post_sigma_c
    )
  },
  numeric(n_post_plot)
)

posterior_f_rep_all <- vapply(
  1:14,
  function(j) {
    10^rnorm(
      n_post_plot,
      mean = log10(
        pmax(post_f_hat[, j], 1e-12)
      ),
      sd = post_sigma_f
    )
  },
  numeric(n_post_plot)
)

colnames(posterior_y_rep) <- paste0(
  "y_rep[",
  seq_len(data_list$NY),
  "]"
)

colnames(posterior_f_rep_all) <- paste0(
  "f_rep[",
  1:14,
  "]"
)


dim(prior_y_rep)
dim(prior_f_rep_all)

dim(posterior_y_rep)
dim(posterior_f_rep_all)





# ============================================================
# 7. Summarise predictive distributions
# ============================================================

# Extract the latest observation without assuming whether time is stored in
# rows or columns.  This is the assumption that caused the previous f_obs
# length error.
extract_latest_vector <- function(
    x,
    accepted_sizes,
    n_time = NULL,
    object_name = "object"
) {
  
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  
  dims <- dim(x)
  
  if (is.null(dims)) {
    
    values <- as.numeric(x)
    
    if (!(length(values) %in% accepted_sizes)) {
      stop(
        object_name,
        " is a vector of length ",
        length(values),
        "; expected one of: ",
        paste(accepted_sizes, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    
    return(values)
  }
  
  if (length(dims) != 2L) {
    stop(
      object_name,
      " must be a vector or a two-dimensional matrix/data frame; dimensions are ",
      paste(dims, collapse = " x "),
      ".",
      call. = FALSE
    )
  }
  
  nr <- dims[1]
  nc <- dims[2]
  
  # Prefer the orientation that agrees with N_t when it is informative.
  if (
    !is.null(n_time) &&
    nr == n_time &&
    nc %in% accepted_sizes
  ) {
    return(as.numeric(x[nr, , drop = TRUE]))
  }
  
  if (
    !is.null(n_time) &&
    nc == n_time &&
    nr %in% accepted_sizes
  ) {
    return(as.numeric(x[, nc, drop = TRUE]))
  }
  
  # Fallback when there is only one plausible variable dimension.
  if (nc %in% accepted_sizes && !(nr %in% accepted_sizes)) {
    return(as.numeric(x[nr, , drop = TRUE]))
  }
  
  if (nr %in% accepted_sizes && !(nc %in% accepted_sizes)) {
    return(as.numeric(x[, nc, drop = TRUE]))
  }
  
  # A one-row or one-column object is also unambiguous.
  if (nr == 1L && nc %in% accepted_sizes) {
    return(as.numeric(x[1, , drop = TRUE]))
  }
  
  if (nc == 1L && nr %in% accepted_sizes) {
    return(as.numeric(x[, 1, drop = TRUE]))
  }
  
  stop(
    "Cannot determine the orientation of ",
    object_name,
    " with dimensions ",
    paste(dims, collapse = " x "),
    ". Expected the variable dimension to have size ",
    paste(accepted_sizes, collapse = " or "),
    ".",
    call. = FALSE
  )
}


# Observed concentrations at the final observation time
observed_y <- extract_latest_vector(
  data_list$y_obs,
  accepted_sizes = data_list$NY,
  n_time = data_list$N_t,
  object_name = "data_list$y_obs"
)

stopifnot(
  length(observed_y) == data_list$NY,
  length(observed_y) == ncol(prior_y_rep),
  length(observed_y) == ncol(posterior_y_rep)
)



summarise_predictive <- function(
    draws,
    observed,
    labels,
    distribution
) {
  
  relative_draws <- sweep(
    draws,
    MARGIN = 2,
    STATS = observed,
    FUN = "/"
  )
  
  quantiles <- t(
    apply(
      relative_draws,
      MARGIN = 2,
      quantile,
      probs = c(
        0.05,
        0.25,
        0.5,
        0.75,
        0.95
      ),
      na.rm = TRUE
    )
  )
  
  data.frame(
    Variable = labels,
    Distribution = distribution,
    Lower90 = quantiles[, 1],
    Lower50 = quantiles[, 2],
    Median = quantiles[, 3],
    Upper50 = quantiles[, 4],
    Upper90 = quantiles[, 5]
  )
}


if (data_list$NY == 13) {
  
  concentration_labels <- c(
    "Pyruvate",
    "Oxaloacetate",
    "Acetyl-CoA",
    "Citrate",
    "cis-Aconitate",
    "Isocitrate",
    "\u03b1-Ketoglutarate",
    "Succinyl-CoA",
    "Succinate",
    "Fumarate",
    "Acetyl-phosphate",
    "NAD\u207A",
    "NADH"
  )
  
} else {
  
  concentration_labels <- paste0(
    "Concentration ",
    seq_len(data_list$NY)
  )
}

flux_labels <- c(
  "vCS",
  "vACO1",
  "vACO2",
  "vIDH",
  "vAKGDH",
  "vSCS",
  "vSDH",
  "vFUM",
  "vMDH",
  "vPTA",
  "vACKA"
)



concentration_summary <- rbind(
  summarise_predictive(
    prior_y_rep,
    observed_y,
    concentration_labels,
    "Prior predictive"
  ),
  summarise_predictive(
    posterior_y_rep,
    observed_y,
    concentration_labels,
    "Posterior predictive"
  )
)

# Extract the 11 observed model fluxes:
# CS, ACO1, ACO2, IDH, AKGDH, SCS, SDH, FUM, MDH, PTA, ACKA

prior_f_rep <- prior_f_rep_all[
  ,
  3:13,
  drop = FALSE
]

# Posterior f_rep may contain either the 11 observed fluxes, the 13 reaction
# fluxes (PK through ACKA), or all 14 outputs including the redox flux.
if (ncol(posterior_f_rep_all) == 11) {
  
  posterior_f_rep <- posterior_f_rep_all
  
} else if (ncol(posterior_f_rep_all) %in% c(13, 14)) {
  
  posterior_f_rep <- posterior_f_rep_all[
    ,
    3:13,
    drop = FALSE
  ]
  
} else {
  
  stop(
    "Unexpected number of columns in posterior_f_rep_all: ",
    ncol(posterior_f_rep_all)
  )
}

# Extract the latest observed flux vector.  Depending on how data_list was
observed_f_tca <- as.numeric(
  data_list$f_obs_tca
)

observed_f_ferm <- as.numeric(
  data_list$f_obs_ferm
)

stopifnot(
  length(observed_f_tca) == 1L,
  length(observed_f_ferm) == 1L
)

observed_f <- c(
  rep(observed_f_tca, 9),
  rep(observed_f_ferm, 2)
)

names(observed_f) <- flux_labels

stopifnot(
  ncol(prior_f_rep) == 11L,
  ncol(posterior_f_rep) == 11L,
  length(observed_f) == 11L
)



flux_summary <- rbind(
  summarise_predictive(
    prior_f_rep,
    observed_f,
    flux_labels,
    "Prior predictive"
  ),
  summarise_predictive(
    posterior_f_rep,
    observed_f,
    flux_labels,
    "Posterior predictive"
  )
)

concentration_summary$Variable <- factor(
  concentration_summary$Variable,
  levels = rev(concentration_labels)
)

flux_summary$Variable <- factor(
  flux_summary$Variable,
  levels = rev(flux_labels)
)



library(ggplot2)
library(dplyr)

prior_fill <- "grey65"
posterior_fill <- "#56B4E9"
ci_colour <- "#005A9C"


# ============================================================
# Convert a predictive matrix to long format
# ============================================================

predictive_matrix_to_long <- function(
    draws,
    labels,
    distribution
) {
  
  # Drop the posterior::draws_matrix class completely.  Otherwise the
  # value column can retain that class and fail when rows are combined.
  draws <- matrix(
    as.numeric(draws),
    nrow = nrow(draws),
    ncol = ncol(draws),
    dimnames = dimnames(draws)
  )
  
  stopifnot(
    ncol(draws) == length(labels)
  )
  
  data.frame(
    Variable = factor(
      rep(
        labels,
        each = nrow(draws)
      ),
      levels = labels
    ),
    Distribution = distribution,
    value = as.numeric(draws)
  )
}


# ============================================================
# Predictive-density plotting function
# ============================================================

make_predictive_density_plot <- function(
    prior_draws,
    posterior_draws,
    observed,
    labels,
    x_label,
    ncol = 5,
    x_limits = NULL,
    facet_scales = "free_x",
    interval_prob = 0.90,
    display_prob = 0.995
) {
  
  prior_draws <- matrix(
    as.numeric(prior_draws),
    nrow = nrow(prior_draws),
    ncol = ncol(prior_draws),
    dimnames = dimnames(prior_draws)
  )
  
  posterior_draws <- matrix(
    as.numeric(posterior_draws),
    nrow = nrow(posterior_draws),
    ncol = ncol(posterior_draws),
    dimnames = dimnames(posterior_draws)
  )
  observed <- as.numeric(observed)
  
  stopifnot(
    ncol(prior_draws) == length(labels),
    ncol(posterior_draws) == length(labels),
    length(observed) == length(labels),
    all(observed > 0)
  )
  
  
  # ----------------------------------------------------------
  # Long-format predictive draws
  # ----------------------------------------------------------
  
  df_all <- rbind(
    predictive_matrix_to_long(
      prior_draws,
      labels,
      "Prior predictive"
    ),
    predictive_matrix_to_long(
      posterior_draws,
      labels,
      "Posterior predictive"
    )
  )
  
  # Remove non-positive/non-finite draws before log transformation
  df_all <- df_all %>%
    filter(
      is.finite(value),
      value > 1e-10
    ) %>%
    mutate(
      log10_value = log10(value)
    )
  
  # Use only the central display_prob of the combined prior/posterior draws
  # to determine each panel's x range.  This prevents a few extreme draws
  # from creating large empty regions; it does not change the 90% PI below.
  display_tail <- (1 - display_prob) / 2
  
  display_bounds <- df_all %>%
    group_by(Variable) %>%
    summarise(
      display_lower = quantile(
        log10_value,
        probs = display_tail,
        na.rm = TRUE
      ),
      display_upper = quantile(
        log10_value,
        probs = 1 - display_tail,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  df_all <- df_all %>%
    left_join(
      display_bounds,
      by = "Variable"
    ) %>%
    filter(
      log10_value >= display_lower,
      log10_value <= display_upper
    )
  
  
  # ----------------------------------------------------------
  # Observed values
  # ----------------------------------------------------------
  
  df_observed <- data.frame(
    Variable = factor(
      labels,
      levels = labels
    ),
    observed_value = observed,
    observed_log10 = log10(observed)
  )
  
  
  # ----------------------------------------------------------
  # Posterior predictive median and central PI
  # ----------------------------------------------------------
  
  posterior_for_pi <- posterior_draws
  posterior_for_pi[
    !is.finite(posterior_for_pi) |
      posterior_for_pi <= 1e-10
  ] <- NA_real_
  
  posterior_log10 <- log10(posterior_for_pi)
  
  interval_tail <- (1 - interval_prob) / 2
  
  posterior_quantiles <- t(
    apply(
      posterior_log10,
      MARGIN = 2,
      quantile,
      probs = c(
        interval_tail,
        0.5,
        1 - interval_tail
      ),
      na.rm = TRUE
    )
  )
  
  df_pi <- data.frame(
    Variable = factor(
      labels,
      levels = labels
    ),
    pi_lower = posterior_quantiles[, 1],
    median_value = posterior_quantiles[, 2],
    pi_upper = posterior_quantiles[, 3]
  )
  
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  pi_legend_label <- paste0(
    "Posterior ",
    round(100 * interval_prob),
    "% PI"
  )
  
  ggplot() +
    
    # Light-blue posterior predictive interval band
    geom_rect(
      data = df_pi,
      mapping = aes(
        xmin = pi_lower,
        xmax = pi_upper,
        ymin = -Inf,
        ymax = Inf,
        fill = pi_legend_label
      ),
      alpha = 0.32,
      colour = NA,
      inherit.aes = FALSE
    ) +
    
    # Prior predictive density
    geom_density(
      data = df_all %>%
        filter(
          Distribution == "Prior predictive"
        ),
      mapping = aes(
        x = log10_value,
        y = after_stat(density),
        fill = Distribution
      ),
      colour = prior_fill,
      linewidth = 0.35,
      alpha = 0.30,
      adjust = 1,
      na.rm = TRUE
    ) +
    
    # Posterior predictive density
    geom_density(
      data = df_all %>%
        filter(
          Distribution == "Posterior predictive"
        ),
      mapping = aes(
        x = log10_value,
        y = after_stat(density),
        fill = Distribution
      ),
      colour = posterior_fill,
      linewidth = 0.45,
      alpha = 0.45,
      adjust = 1,
      na.rm = TRUE
    ) +
    
    # Posterior predictive median: blue dashed vertical line
    geom_vline(
      data = df_pi,
      mapping = aes(
        xintercept = median_value,
        colour = "Posterior median",
        linetype = "Posterior median"
      ),
      linewidth = 0.65
    ) +
    
    # Observed value: orange solid vertical line
    geom_vline(
      data = df_observed,
      mapping = aes(
        xintercept = observed_log10,
        colour = "Observed",
        linetype = "Observed"
      ),
      linewidth = 0.65
    ) +
    
    facet_wrap(
      facets = vars(Variable),
      ncol = ncol,
      scales = facet_scales
    ) +
    
    scale_fill_manual(
      name = NULL,
      values = c(
        "Prior predictive" = prior_fill,
        "Posterior predictive" = posterior_fill,
        "Posterior 90% PI" = "#D9EEF8"
      )
    ) +
    
    scale_colour_manual(
      name = NULL,
      values = c(
        "Posterior median" = ci_colour,
        "Observed" = "#E66101"
      )
    ) +
    
    scale_linetype_manual(
      name = NULL,
      values = c(
        "Posterior median" = "dashed",
        "Observed" = "solid"
      )
    ) +
    
    scale_y_continuous(
      breaks = scales::breaks_pretty(
        n = 3
      ),
      expand = expansion(
        mult = c(0.03, 0.08)
      )
    ) +
    
    coord_cartesian(
      xlim = x_limits,
      clip = "on"
    ) +
    
    labs(
      x = x_label,
      y = "Predictive density"
    ) +
    
    guides(
      fill = guide_legend(
        order = 1,
        override.aes = list(
          alpha = c(0.30, 0.45, 0.32)
        )
      ),
      colour = guide_legend(
        order = 2,
        override.aes = list(linewidth = 0.8)
      ),
      linetype = guide_legend(
        order = 2
      )
    ) +
    
    theme_classic(
      base_size = 9,
      base_family = "sans"
    ) +
    
    theme(
      legend.position = "top",
      legend.justification = "center",
      
      strip.background = element_blank(),
      
      strip.text = element_text(
        size = 9,
        face = "plain",
        colour = "black"
      ),
      
      axis.title = element_text(
        size = 9,
        colour = "black"
      ),
      
      axis.text = element_text(
        size = 7,
        colour = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.35,
        colour = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.4,
        colour = "black"
      ),
      
      panel.spacing = grid::unit(
        1,
        "lines"
      ),
      
      plot.margin = margin(
        t = 5,
        r = 5,
        b = 5,
        l = 5
      )
    )
}



concentration_density_plot <- make_predictive_density_plot(
  prior_draws = prior_y_rep,
  posterior_draws = posterior_y_rep,
  observed = observed_y,
  labels = concentration_labels,
  x_label = expression(
    log[10] * " predicted concentration (mM)"
  ),
  ncol = 5,
  x_limits = NULL,
  facet_scales = "free_x",
  interval_prob = 0.90,
  display_prob = 0.995
)



flux_density_plot <- make_predictive_density_plot(
  prior_draws = prior_f_rep,
  posterior_draws = posterior_f_rep,
  observed = observed_f,
  labels = flux_labels,
  x_label = expression(
    log[10] * " flux (mM " * s^{-1} * ")"
  ),
  ncol = 5,
  x_limits = NULL,
  facet_scales = "free_x",
  interval_prob = 0.90,
  display_prob = 0.995
)

print(concentration_density_plot)
print(flux_density_plot)


library(patchwork)

final_fig <- (
  concentration_density_plot /
    flux_density_plot
) +
  plot_layout(
    heights = c(1, 1),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "(",
    tag_suffix = ")"
  )

final_fig <- final_fig &
  theme(
    legend.position = "top",
    
    plot.tag = element_text(
      size = 10,
      face = "bold",
      colour = "black"
    ),
    
    plot.tag.position = c(0, 1)
  )

print(final_fig)


save_journal_plot(concentration_density_plot, "tcaferm_real_priorposterior_met_ppc0921", width = 20, height = 10)
save_journal_plot(flux_density_plot, "tcaferm_real_priorposterior_flux_ppc0921", width = 20, height = 10)
save_journal_plot(final_fig, "tcaferm_real_priorposterior_final_fig0921", width = 20, height = 20)
