model {
  // No parameters and no likelihood: use fixed_param sampling.
}

generated quantities {
  real pdh_raw     = normal_rng(-2, 1);
  real aco_raw     = normal_rng(-2, 1);
  real akd_raw     = normal_rng(-2, 1);
  real scs_raw     = normal_rng(log10(0.012340987), 0.5);
  real sd_raw      = normal_rng(-2, 1);
  real pk_raw      = normal_rng(-2, 1);
  real cs_raw      = normal_rng(-2, 1);
  real id_raw      = normal_rng(-2, 1);
  real pta_raw     = normal_rng(log10(0.000887043), 0.5);

  // Retained for compatibility with the fitted model's parameter block.
  // It does not affect predictions while theta[10] is fixed at 5.
  real acetate_raw = normal_rng(-1, 1);

  real md_raw      = normal_rng(log10(0.023603447), 0.5);
  real f_raw       = normal_rng(log10(0.001278269), 0.5);
  real ack_raw     = normal_rng(-2, 1);
  real k_raw       = normal_rng(log10(83), 0.5);
  real alpha_raw   = normal_rng(log10(0.0545 / 2.8), 0.5);

  real sigma_c_raw = normal_rng(log10(0.2), 0.15);
  real sigma_f_raw = normal_rng(log10(0.2), 0.15);
}
