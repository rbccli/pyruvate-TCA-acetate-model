library(dplyr)
library(tidyr)
library(ggplot2)
library(posterior)
library(scales)

y_obs <- c(
  pyruvate      = 2.439091e+01,
  oxaloacetate  = 6.901385e-03,
  acetyl_coa    = 1.430352e-02,
  citrate       = 5.235775e-02,
  cis_aconitate = 1.112080e-03,
  isocitrate    = 3.264758e-05,
  aketo         = 3.588200e-05,
  succinyl_coa  = 1.183278e-04,
  succinate     = 9.684217e-06,
  fumarate      = 1.072823e+00,
  acetyl_p      = 1.153468e-03,
  nad           = 3.169765e+00,
  nadh          = 6.335408e-02
)

f_obs_ferm <- 0.0609363 
f_obs_tca <- 0.001873996


met_labels <- c(
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

flux_labels <- c(
  "vPK", "vPDH", "vCS", "vACO1", "vACO2",
  "vIDH", "vAKGDH", "vSCS", "vSDH", "vFUM",
  "vMDH", "vPTA", "vACK"
)


draws <- tcaferm_fit_simulation$draws()


# ============================================================
# Final journal theme
# ============================================================

theme_ppc_final <- function(base_size = 8.5,
                            base_family = "sans") {
  
  theme_classic(
    base_size = base_size,
    base_family = base_family
  ) +
    theme(
      axis.title = element_text(
        size = 9,
        colour = "black"
      ),
      
      axis.text = element_text(
        size = 7.5,
        colour = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.35,
        colour = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.30,
        colour = "black"
      ),
      
      axis.ticks.length = unit(1.5, "mm"),
      
      strip.background = element_blank(),
      
      strip.text = element_text(
        size = 8.5,
        face = "plain",
        colour = "black",
        margin = margin(b = 2.5)
      ),
      
      panel.spacing.x = unit(0.48, "cm"),
      panel.spacing.y = unit(0.40, "cm"),
      
      plot.margin = margin(
        t = 4,
        r = 5,
        b = 4,
        l = 5
      )
    )
}


# ============================================================
# Select a visually sensible common tick interval
# ============================================================

nice_tick_step <- function(target_step) {
  
  magnitude <- 10^floor(log10(target_step))
  
  candidates <- magnitude * c(
    1, 2, 2.5, 5, 10
  )
  
  candidates[
    which.min(abs(candidates - target_step))
  ]
}


# ============================================================
# Reusable PPC histogram function
# ============================================================

make_ppc_histogram <- function(
    draws,
    prefix,
    panel_labels,
    observed,
    x_title,
    ncol = 5,
    interval_prob = 0.9,
    display_prob = 0.955,
    number_of_bins = 28,
    posterior_colour = "#9ECAE1",
    interval_colour = "#4292C6",
    observed_colour = "#D55E00") {
  
  n_panels <- length(panel_labels)
  
  variable_names <- sprintf(
    "%s[%d]",
    prefix,
    seq_len(n_panels)
  )
  
  draws_df <- posterior::as_draws_df(draws)
  
  missing_variables <- setdiff(
    variable_names,
    names(draws_df)
  )
  
  if (length(missing_variables) > 0) {
    stop(
      "The following variables are missing from the posterior draws: ",
      paste(missing_variables, collapse = ", ")
    )
  }
  
  # ----------------------------------------------------------
  # Posterior draws in long format
  # ----------------------------------------------------------
  
  posterior_long <- draws_df |>
    dplyr::select(all_of(variable_names)) |>
    mutate(draw = row_number()) |>
    pivot_longer(
      cols = all_of(variable_names),
      names_to = "variable",
      values_to = "value"
    ) |>
    mutate(
      panel_id = match(variable, variable_names),
      
      panel = factor(
        panel_labels[panel_id],
        levels = panel_labels
      )
    ) |>
    filter(
      is.finite(value),
      value > 0
    ) |>
    mutate(
      log_value = log10(value)
    )
  
  # ----------------------------------------------------------
  # Observed values
  # ----------------------------------------------------------
  
  observed_df <- tibble(
    panel_id = seq_len(n_panels),
    observed = as.numeric(observed)
  ) |>
    mutate(
      panel = factor(
        panel_labels[panel_id],
        levels = panel_labels
      ),
      
      log_observed = log10(observed)
    ) |>
    filter(
      is.finite(observed),
      observed > 0
    )
  
  # ----------------------------------------------------------
  # Central 95% posterior predictive intervals
  # ----------------------------------------------------------
  
  interval_tail <- (1 - interval_prob) / 2
  
  interval_df <- posterior_long |>
    group_by(panel) |>
    summarise(
      interval_low = quantile(
        log_value,
        interval_tail,
        na.rm = TRUE
      ),
      
      posterior_median = median(
        log_value,
        na.rm = TRUE
      ),
      
      interval_high = quantile(
        log_value,
        1 - interval_tail,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    )
  
  # ----------------------------------------------------------
  # Display range
  #
  # Use the central 99.8% of the draws to prevent a small
  # number of extreme draws from setting the complete x-range.
  # The observed value is always included.
  # ----------------------------------------------------------
  
  display_tail <- (1 - display_prob) / 2
  
  range_df <- posterior_long |>
    group_by(panel) |>
    summarise(
      posterior_low = quantile(
        log_value,
        display_tail,
        na.rm = TRUE
      ),
      
      posterior_high = quantile(
        log_value,
        1 - display_tail,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) |>
    left_join(
      observed_df |>
        select(panel, log_observed),
      by = "panel"
    ) |>
    mutate(
      range_low = pmin(
        posterior_low,
        log_observed
      ),
      
      range_high = pmax(
        posterior_high,
        log_observed
      ),
      
      panel_midpoint = (
        range_low + range_high
      ) / 2,
      
      required_width = (
        range_high - range_low
      )
    )
  
  # Same log10 width for every panel
  common_width <- max(
    range_df$required_width,
    na.rm = TRUE
  ) * 1.10
  
  if (!is.finite(common_width) || common_width <= 0) {
    stop("Unable to calculate a valid plotting range.")
  }
  
  # Same major-tick spacing for every panel
  tick_step <- nice_tick_step(
    common_width / 4
  )
  
  x_break_function <- function(limits) {
    
    first_break <- ceiling(
      limits[1] / tick_step
    ) * tick_step
    
    last_break <- floor(
      limits[2] / tick_step
    ) * tick_step
    
    if (first_break > last_break) {
      return(numeric(0))
    }
    
    seq(
      first_break,
      last_break,
      by = tick_step
    )
  }
  
  bounds_df <- range_df |>
    mutate(
      xmin = panel_midpoint - common_width / 2,
      xmax = panel_midpoint + common_width / 2
    ) |>
    select(
      panel,
      xmin,
      xmax
    )
  
  # Clip only the extreme plotting tails
  posterior_clipped <- posterior_long |>
    left_join(
      bounds_df,
      by = "panel"
    ) |>
    filter(
      log_value >= xmin,
      log_value <= xmax
    )
  
  # Invisible anchors force every facet to have the same width
  xlim_df <- bounds_df |>
    pivot_longer(
      cols = c(xmin, xmax),
      names_to = "bound",
      values_to = "log_value"
    ) |>
    mutate(y = 0)
  
  histogram_binwidth <- (
    common_width / number_of_bins
  )
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  ggplot(
    posterior_clipped,
    aes(x = log_value)
  ) +
    
    geom_blank(
      data = xlim_df,
      aes(
        x = log_value,
        y = y
      ),
      inherit.aes = FALSE
    ) +
    
    geom_rect(
      data = interval_df,
      aes(
        xmin = interval_low,
        xmax = interval_high,
        ymin = 0,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = interval_colour,
      alpha = 0.11
    ) +
    
    geom_histogram(
      aes(
        y = after_stat(density)
      ),
      binwidth = histogram_binwidth,
      boundary = 0,
      fill = posterior_colour,
      alpha = 0.80,
      colour = NA
    ) +
    
    geom_vline(
      data = observed_df,
      aes(
        xintercept = log_observed
      ),
      linewidth = 0.65,
      colour = observed_colour
    ) +
    
    facet_wrap(
      vars(panel),
      scales = "free_x",
      ncol = ncol,
      drop = FALSE
    ) +
    
    scale_x_continuous(
      breaks = x_break_function,
      labels = label_number(
        accuracy = tick_step,
        trim = TRUE
      ),
      expand = expansion(mult = 0)
    ) +
    
    scale_y_continuous(
      limits = c(0, NA),
      breaks = breaks_extended(n = 3),
      labels = label_number(
        accuracy = 0.1,
        trim = TRUE
      ),
      expand = expansion(
        mult = c(0, 0.08)
      )
    ) +
    
    labs(
      x = x_title,
      y = "Posterior predictive density"
    ) +
    
    theme_ppc_final()
}


# ============================================================
# (a) Metabolite PPC
# ============================================================

p_met_ppc <- make_ppc_histogram(
  draws = draws,
  prefix = "y_rep",
  panel_labels = met_labels,
  observed = y_obs,
  
  x_title = expression(
    log[10] * " concentration (mM)"
  ),
  
  ncol = 5
)


# ============================================================
# (b) Flux PPC
# ============================================================

f_obs <- c(
  f_obs_ferm + f_obs_tca,  # vPK
  f_obs_ferm + f_obs_tca,  # vPDH
  
  rep(
    f_obs_tca,
    9
  ),                       # vCS to vMDH
  
  f_obs_ferm,              # vPTA
  f_obs_ferm               # vACK
)

stopifnot(
  length(f_obs) == length(flux_labels)
)

p_flux_ppc <- make_ppc_histogram(
  draws = draws,
  prefix = "f_rep",
  panel_labels = flux_labels,
  observed = f_obs,
  
  x_title = expression(
    log[10] * " flux (mM " * s^{-1} * ")"
  ),
  
  ncol = 5
)


# ============================================================
# Combined figure
# ============================================================

final_fig <- (
  p_met_ppc /
    p_flux_ppc
) +
  plot_layout(
    heights = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "(",
    tag_suffix = ")"
  )

final_fig <- final_fig &
  theme(
    plot.tag = element_text(
      size = 10,
      colour = "black"
    ),
    
    plot.tag.position = c(0, 1)
  )

final_fig



save_journal_plot(p_met_ppc, "tcaferm_simulation_met_ppc", width = 20, height = 10)
save_journal_plot(p_flux_ppc, "tcaferm_simulation_flux_ppc", width = 20, height = 10)
save_journal_plot(final_fig, "tcaferm_simulation_final_fig", width = 20, height = 15)





