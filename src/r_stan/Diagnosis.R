# ============================================================
# MCMC convergence and sampling diagnostics
# ============================================================

library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)


# ============================================================
# 1. Parameters sampled in Stan
# ============================================================

raw_params <- c(
  "pdh_raw",
  "akd_raw",
  "scs_raw",
  "aco_raw",
  "sd_raw",
  "pk_raw",
  "cs_raw",
  "id_raw",
  "pta_raw",
  "md_raw",
  "f_raw",
  "ack_raw",
  "alpha_raw",
  "k_raw",
  "sigma_c",
  "sigma_f"
)


# ============================================================
# 2. Check that all parameters exist
# ============================================================

available_variables <- posterior::variables(
  tca_ferm_fit_real$draws()
)

missing_params <- setdiff(
  raw_params,
  available_variables
)

if (length(missing_params) > 0) {
  stop(
    paste(
      "The following parameters were not found:",
      paste(missing_params, collapse = ", ")
    )
  )
}

message("All requested parameters were found.")


# ============================================================
# 3. Extract post-warm-up draws
# ============================================================

# By default, fit$draws() returns post-warm-up draws only.

trace_draws <- tca_ferm_fit_real$draws(
  variables = raw_params,
  format = "draws_array"
)


# ============================================================
# 4. Split parameters into two pages
# ============================================================

parameter_groups <- list(
  page_1 = raw_params[1:7],
  page_2 = raw_params[8:14]
)


# ============================================================
# 5. Create trace plots
# ============================================================

bayesplot::color_scheme_set("mix-blue-red")

trace_plots <- lapply(
  seq_along(parameter_groups),
  function(i) {
    
    current_parameters <- parameter_groups[[i]]
    
    bayesplot::mcmc_trace(
      trace_draws,
      pars = current_parameters,
      facet_args = list(
        ncol = 2,
        scales = "free_y"
      ),
      size = 0.35,
      alpha = 0.8
    ) +
      labs(
        title = paste0(
          "Post-warm-up MCMC trace plots (",
          i,
          "/",
          length(parameter_groups),
          ")"
        ),
        x = "Post-warm-up iteration",
        y = NULL,
        colour = "Chain"
      ) +
      theme_bw(base_size = 11) +
      theme(
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(
          size = 10,
          face = "bold"
        ),
        panel.grid.minor = element_blank(),
        plot.title = element_text(
          size = 12,
          face = "bold"
        )
      )
  }
)


# Display plots in RStudio

for (p in trace_plots) {
  print(p)
}


# ============================================================
# 6. Save trace plots as a multi-page PDF
# ============================================================

trace_file <- "mcmc_traceplots.pdf"

grDevices::cairo_pdf(
  filename = trace_file,
  width = 10,
  height = 10,
  onefile = TRUE
)

for (p in trace_plots) {
  print(p)
}

grDevices::dev.off()

message(
  "Trace plots saved to: ",
  normalizePath(trace_file)
)


# ============================================================
# 7. Produce parameter-level diagnostic table
# ============================================================

diagnostic_summary <- tca_ferm_fit_real$summary(
  variables = raw_params
)

diagnostic_table <- diagnostic_summary[
  ,
  c(
    "variable",
    "rhat",
    "ess_bulk",
    "ess_tail"
  )
]

names(diagnostic_table) <- c(
  "Parameter",
  "Rhat",
  "Bulk_ESS",
  "Tail_ESS"
)

# Round values for display

diagnostic_table$Rhat <- formatC(
  diagnostic_table$Rhat,
  format = "f",
  digits = 3
)

diagnostic_table$Bulk_ESS <- formatC(
  diagnostic_table$Bulk_ESS,
  format = "f",
  digits = 3
)

diagnostic_table$Tail_ESS <- formatC(
  diagnostic_table$Tail_ESS,
  format = "f",
  digits = 3
)

print(
  diagnostic_table,
  row.names = FALSE
)


# Save diagnostic table

write.csv(
  diagnostic_table,
  file = "mcmc_parameter_diagnostics.csv",
  row.names = FALSE
)

message(
  "Diagnostic table saved to: ",
  normalizePath("mcmc_parameter_diagnostics.csv")
)


# ============================================================
# 8. Summarise R-hat and ESS
# ============================================================

cat("\nParameter-level diagnostic summary\n")
cat("----------------------------------\n")

cat(
  "Maximum R-hat:",
  max(
    diagnostic_summary$rhat,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Minimum bulk ESS:",
  min(
    diagnostic_summary$ess_bulk,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Minimum tail ESS:",
  min(
    diagnostic_summary$ess_tail,
    na.rm = TRUE
  ),
  "\n"
)

# Identify parameters with the least efficient sampling

lowest_bulk_parameter <- diagnostic_summary$variable[
  which.min(diagnostic_summary$ess_bulk)
]

lowest_tail_parameter <- diagnostic_summary$variable[
  which.min(diagnostic_summary$ess_tail)
]

cat(
  "Parameter with lowest bulk ESS:",
  lowest_bulk_parameter,
  "\n"
)

cat(
  "Parameter with lowest tail ESS:",
  lowest_tail_parameter,
  "\n"
)


# ============================================================
# 9. Divergences and maximum tree depth
# ============================================================

sampler_diagnostics <- tca_ferm_fit_real$sampler_diagnostics(
  format = "draws_matrix"
)

sampler_diagnostics <- as.matrix(
  sampler_diagnostics
)

number_of_postwarmup_draws <- nrow(
  sampler_diagnostics
)

number_of_divergences <- sum(
  sampler_diagnostics[, "divergent__"]
)

divergence_percentage <- (
  100 *
    number_of_divergences /
    number_of_postwarmup_draws
)

# Change this if a different max_treedepth was used.

max_treedepth_used <- 12

number_at_max_treedepth <- sum(
  sampler_diagnostics[, "treedepth__"] >=
    max_treedepth_used
)

cat("\nSampler-level diagnostic summary\n")
cat("--------------------------------\n")

cat(
  "Post-warm-up draws:",
  number_of_postwarmup_draws,
  "\n"
)

cat(
  "Divergent transitions:",
  number_of_divergences,
  sprintf(
    "(%.3f%%)\n",
    divergence_percentage
  )
)

cat(
  "Transitions reaching maximum tree depth:",
  number_at_max_treedepth,
  "\n"
)


# ============================================================
# 10. CmdStanR diagnostic summary
# ============================================================

tca_ferm_fit_real$diagnostic_summary()
