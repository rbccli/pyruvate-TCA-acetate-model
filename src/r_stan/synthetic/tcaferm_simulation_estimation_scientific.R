# install.packages("deSolve")
library(deSolve)


# inhibition functions
atp_inhibition <- function(atp_con, ki_atp) {
  return(1 / (1 + atp_con / ki_atp))
  
}

nadh_inhibition <- function(nadh_con, ki_nadh) {
  return(1 / (1 + nadh_con / ki_nadh))
}

# ACO (aconitase)
compute_v_aco <- function(citrate_con, cis_aconitate_con, isocitrate_con,
                          km_cit, km_cis1, km_cis2, km_iso,
                          enzyme_aco, kf1, kr1, kf2, kr2) {
  cit_tilde <- citrate_con / km_cit
  cis_tilde1 <- cis_aconitate_con / km_cis1
  cis_tilde2 <- cis_aconitate_con / km_cis2
  iso_tilde <- isocitrate_con / km_iso
  
  denom1 <- 1 + cit_tilde + cis_tilde1
  denom2 <- 1 + cis_tilde2 + iso_tilde
  
  v1 <- (kf1 * cit_tilde - kr1 * cis_tilde1) / denom1 * enzyme_aco
  v2 <- (kf2 * cis_tilde2 - kr2 * iso_tilde) / denom2 * enzyme_aco
  
  return(c(v1, v2))
}


# AKD (alpha-ketoglutarate dehydrogenase)
compute_v_akd <- function(a_keto_con, nad_con, coa_sh_con, succinyl_coa_con,
                          nadh_con, co2_con,
                          km_co2, km_ak, km_nad, km_nadh,
                          km_coa, km_succinyl_coa,
                          enzyme_akd, kf_akd, kr_akd) {
  ak_tilde  <- a_keto_con / km_ak
  nad_tilde <- nad_con / km_nad
  coa_tilde <- coa_sh_con / km_coa
  succ_tilde <- succinyl_coa_con / km_succinyl_coa
  nadh_tilde <- nadh_con / km_nadh
  co2_tilde <- co2_con / km_co2
  
  denom <- (1 + nad_tilde) * (1 + ak_tilde) * (1 + coa_tilde) +
    (1 + succ_tilde) * (1 + nadh_tilde) * (1 + co2_tilde) - 1
  
  v <- (kf_akd * ak_tilde * coa_tilde * nad_tilde - kr_akd * succ_tilde * nadh_tilde * co2_tilde) / denom * enzyme_akd
  
  return(v)
}


# CS (citrate synthase)
compute_v_cs <- function(oxaloacetate_con, acetyl_coa_con, citrate_con, coa_sh_con,
                         nadh_con,
                         km_ox, km_acc, km_cit, km_coa_sh,
                         enzyme_cs, kf_cs, kr_cs, ki_nadh) {
  ox_tilde  <- oxaloacetate_con / km_ox
  acc_tilde <- acetyl_coa_con / km_acc
  cit_tilde <- citrate_con / km_cit
  coa_tilde <- coa_sh_con / km_coa_sh
  
  denom <- (1 + ox_tilde) * (1 + acc_tilde) +
    (1 + cit_tilde) * (1 + coa_tilde) - 1
  
  fwd <- kf_cs * ox_tilde * acc_tilde
  rev <- kr_cs * cit_tilde * coa_tilde
  
  v <- (fwd - rev) / denom * enzyme_cs
  v <- v * nadh_inhibition(nadh_con, ki_nadh)
  
  return(v)
}


# Fumarase
compute_v_fumarase <- function(fumarate_con, malate_con,
                               km_fum, km_mal,
                               enzyme_fum, kf_fum, kr_fum) {
  fum_tilde <- fumarate_con / km_fum
  mal_tilde <- malate_con / km_mal
  
  denom <- 1 + fum_tilde + mal_tilde
  
  v <- (kf_fum * fum_tilde - kr_fum * mal_tilde) / denom * enzyme_fum
  return(v)
}


# ID (isocitrate dehydrogenase)
compute_v_id <- function(isocitrate_con, nad_con, nadh_con,
                         aketo_con, co2_con,
                         km_co2, km_iso, km_nadh, km_nad, km_aketo,
                         enzyme_id, kf_id, kr_id) {
  iso_tilde   <- isocitrate_con / km_iso
  nad_tilde   <- nad_con / km_nad
  nadh_tilde  <- nadh_con / km_nadh
  aketo_tilde <- aketo_con / km_aketo
  co2_tilde   <- co2_con / km_co2
  
  denom <- (1 + iso_tilde) * (1 + nad_tilde) +
    (1 + nadh_tilde) * (1 + aketo_tilde) * (1 + co2_tilde) - 1
  
  v <- (kf_id * iso_tilde * nad_tilde - kr_id * nadh_tilde * aketo_tilde * co2_tilde) / denom * enzyme_id
  
  return(v)
}


# MD (malate dehydrogenase)
compute_v_md <- function(malate_con, nad_con, oaa_con, nadh_con,
                         km_mal, km_nad, km_oaa, km_nadh,
                         enzyme_md, kf_md, kr_md) {
  mal_tilde  <- malate_con / km_mal
  nad_tilde  <- nad_con / km_nad
  oaa_tilde  <- oaa_con / km_oaa
  nadh_tilde <- nadh_con / km_nadh
  
  denom <- (1 + mal_tilde) * (1 + nad_tilde) +
    (1 + oaa_tilde) * (1 + nadh_tilde) - 1
  
  v <- (kf_md * mal_tilde * nad_tilde - kr_md * oaa_tilde * nadh_tilde) / denom * enzyme_md
  
  #return(v)
  return(v)
}


# PDH (pyruvate dehydrogenase)
compute_v_pdh <- function(pyruvate_con, coa_con, nad_con,
                          acetyl_coa_con, nadh_con, co2_con,
                          km_co2, km_pyruvate, km_coa, km_nad,
                          km_acetyl_coa, km_nadh,
                          enzyme_pdh, kf_pdh, kr_pdh, ki_nadh) {
  pyruvate_tilde   <- pyruvate_con / km_pyruvate
  coa_tilde        <- coa_con / km_coa
  nad_tilde        <- nad_con / km_nad
  acetyl_coa_tilde <- acetyl_coa_con / km_acetyl_coa
  nadh_tilde       <- nadh_con / km_nadh
  co2_tilde        <- co2_con / km_co2
  
  denom <- (1 + pyruvate_tilde) * (1 + coa_tilde) * (1 + nad_tilde) +
    (1 + acetyl_coa_tilde) * (1 + nadh_tilde) * (1 + co2_tilde) - 1
  
  v <- (kf_pdh * pyruvate_tilde * coa_tilde * nad_tilde - kr_pdh * acetyl_coa_tilde * nadh_tilde * co2_tilde) / denom * enzyme_pdh
  v <- v * nadh_inhibition(nadh_con, ki_nadh)
  
  return(v)
}


# PK (pyruvate kinase)
compute_v_pk <- function(km_pyruvate, km_atp, km_pep, km_adp,
                         kf_pk, kr_pk, enzyme_pk,
                         pyruvate_con, atp_con, pep_con, adp_con, ki_atp) {
  pyruvate_tilde <- pyruvate_con / km_pyruvate
  pep_tilde      <- pep_con / km_pep
  atp_tilde      <- atp_con / km_atp
  adp_tilde      <- adp_con / km_adp
  
  denom <- (1 + pyruvate_tilde) * (1 + atp_tilde) +
    (1 + pep_tilde) * (1 + adp_tilde) - 1
  
  v <- (kf_pk * pep_tilde * adp_tilde - kr_pk * pyruvate_tilde * atp_tilde) / denom * enzyme_pk
  
  v <- v * atp_inhibition(atp_con, ki_atp)
  return(v)
}


# SCS (succinyl-CoA synthetase)
compute_v_scs <- function(succinyl_coa_con, adp_con, succinate_con,
                          coa_sh_con, atp_con, pi_con,
                          km_succinyl, km_adp, km_succinate,
                          km_coa, km_atp, km_pi,
                          enzyme_scs, kf_scs, kr_scs) {
  sucCoA_tilde <- succinyl_coa_con / km_succinyl
  adp_tilde    <- adp_con / km_adp
  suc_tilde    <- succinate_con / km_succinate
  coa_tilde    <- coa_sh_con / km_coa
  atp_tilde    <- atp_con / km_atp
  pi_tilde    <- pi_con / km_pi
  
  
  denom <- (1 + sucCoA_tilde) * (1 + adp_tilde) * (1 + pi_tilde) +
    (1 + suc_tilde) * (1 + coa_tilde) * (1 + atp_tilde) - 1
  
  v <- (kf_scs * sucCoA_tilde * adp_tilde * pi_tilde -
          kr_scs * suc_tilde * coa_tilde * atp_tilde) / denom * enzyme_scs
  
  return(v)
}


# SD (succinate dehydrogenase)
compute_v_sd <- function(succinate_con, uq_con, fumarate_con, uqh2_con,
                         km_suc, km_uq, km_fum, km_uqh2,
                         enzyme_sd, kf_sd, kr_sd) {
  suc_tilde  <- succinate_con / km_suc
  uq_tilde   <- uq_con / km_uq
  fum_tilde  <- fumarate_con / km_fum
  uqh2_tilde <- uqh2_con / km_uqh2
  
  denom <- (1 + suc_tilde) * (1 + uq_tilde) +
    (1 + fum_tilde) * (1 + uqh2_tilde) - 1
  
  v <- (kf_sd * suc_tilde * uq_tilde -
          kr_sd * fum_tilde * uqh2_tilde) / denom * enzyme_sd
  
  return(v)

}


compute_v_pta <- function(acetyl_coa_con, pi_con, acetyl_p_con, coa_sh_con, 
                          km_acetyl_coa_pta, km_pi_pta, km_acetyl_p_pta, km_coa_pta, 
                          kf_pta, kr_pta, enzyme_pta) {
  acetyl_coa_tilde = acetyl_coa_con / km_acetyl_coa_pta
  pi_tilde = pi_con / km_pi_pta
  acetyl_p_tilde = acetyl_p_con / km_acetyl_p_pta
  coa_tilde = coa_sh_con / km_coa_pta
  
  
  denom = (1 + acetyl_coa_tilde) * (1 + pi_tilde) + (1 + acetyl_p_tilde) * (1 + coa_tilde) - 1
  
  v = (kf_pta * acetyl_coa_tilde * pi_tilde - 
         kr_pta * acetyl_p_tilde * coa_tilde)/ denom * enzyme_pta
  
  
  return(v)
  
}


compute_v_ack <- function(acetyl_p_con, adp_con, acetate_con, atp_con, 
                          km_acetyl_p_ack, km_adp_ack, km_acetate_ack, km_atp_ack, 
                          kf_ack, kr_ack, enzyme_ack){

  
  acetyl_p_tilde = acetyl_p_con / km_acetyl_p_ack
  adp_tilde = adp_con / km_adp_ack
  
  atp_tilde = atp_con / km_atp_ack
  acetate_tilde = acetate_con / km_acetate_ack
  
  
  denom = (1 + acetyl_p_tilde) * (1 + adp_tilde) + (1 + acetate_tilde) * (1 + atp_tilde) - 1
  
  v = (kf_ack * acetyl_p_tilde * adp_tilde - kr_ack * acetate_tilde * atp_tilde) / denom * enzyme_ack
  
  return(v)

}


compute_v_redox <- function(nadh_con, nad_con, k, alpha){
  
  
  v = k * alpha * nad_con - k * nadh_con
  
  return(v)
  
}


## Parameter List
params_true <- list(
  # inihibition constant
  
  ki_nadh_pdh = 0.001,
  ki_nadh_cs = 0.0028,
  
  ki_atp_pk = 4,
  
  
  # michaelis constant
  # pk
  km_pyruvate_pk = 0.0505,
  km_pep_pk = 0.2666,
  km_atp_pk = 16.1752,
  km_adp_pk = 0.1919,
  
  
  # pdh
  km_pyruvate_pdh = 0.3296, 
  km_coa_pdh = 0.0471, 
  km_nad_pdh = 0.0619, 
  km_nadh_pdh = 0.215,
  km_acetyl_coa_pdh = 0.215,
  km_co2_pdh = 0.215,
  
  # pta
  km_acetyl_coa_pta = 0.0491,
  km_pi_pta = 1.9606,
  km_acetyl_p_pta = 0.6297,
  km_coa_pta = 0.0521,
  
  # ack
  km_acetyl_p_ack = 0.0938,
  km_adp_ack = 0.2799,
  km_acetate_ack = 7.9678,
  km_atp_ack = 0.096,
  
  
  # cs
  km_citrate_cs = 0.0018, 
  km_coa_sh_cs = 44914.6131, # 0.0027, 
  km_oxaloacetate_cs = 0.0119, 
  km_acetyl_coa_cs = 0.0521, 
  
  
  # aconitase
  km_citrate_aco = 0.6464, 
  km_cis_aconitate_aco1 = 0.0095,
  km_cis_aconitate_aco2 = 0.2098,
  km_isocitrate_aco = 0.0872, 
  
  # isocitrate dehydrogenase
  km_isocitrate_id = 0.0109, 
  km_nad_id = 1.5211, 
  km_aketo_id = 0.4744, 
  km_nadh_id = 0.007, # 24.3, 
  km_co2_id = 2.4317,
  
  
  
  # akd checked not much different
  km_nad_akd = 0.0512, 
  km_aketo_akd = 0.0026, 
  km_coa_sh_akd = 0.0512, 
  km_nadh_akd = 0.1953, 
  km_succinyl_coa_akd = 0.1953, 
  km_co2_akd = 0.1953,
  
  
  # scs finished check 
  km_coa_sh_scs = 2.0659, 
  km_succinyl_coa_scs = 0.0316, 
  km_adp_scs = 0.3244, 
  km_succinate_scs = 0.1444, 
  km_atp_scs = 0.0426, 
  km_pi_scs = 0.8953,
  
  
  # sd finished check 
  km_succinate_sd = 0.0021,  
  km_fad_sd = 0.002, # q
  km_fumarate_sd = 0.0233, 
  km_fadh2_sd = 12.5427, # qh2
  
  
  # fumarase
  km_fumarate_fumarase = 0.4222, 
  km_malate_fumarase = 0.6905, 
  
  # malate dehydrogenase dont change same ecoli exp
  km_nad_md = 0.418, 
  km_oxaloacetate_md = 0.0244, 
  km_malate_md = 3.8086, 
  km_nadh_md = 0.0301, 
  
  
  
  
  
  enzyme_pk = 0.1,
  enzyme_pdh = 0.1,
  enzyme_cs = 0.1,
  enzyme_aco = 0.1,
  enzyme_id = 0.1,
  enzyme_akd = 0.1,
  enzyme_scs = 0.1,
  enzyme_sd = 0.1,
  enzyme_fumarase = 0.1,
  enzyme_md = 0.1,
  enzyme_pta = 0.1,
  enzyme_ack = 0.1,
  
  
  # forward and reverse k_cat 
  
  kf_pk = 555.4713,
  kr_pk = 0.4124,
  kf_pdh = 36.1192, 
  kr_pdh = 0,
  kf_cs = 141.5176, 
  kr_cs = 5.8632, 
  kf1_aco = 5.7836, 
  kf2_aco = 4.5447,
  kr1_aco = 3.8424, 
  kr2_aco = 1.0233,
  kf_id = 82.0188, 
  kr_id = 14.7374,
  kf_akd = 2.0345, 
  kr_akd = 0.0011,
  kf_scs = 35.8092, 
  kr_scs = 60.5031, 
  kf_sd = 32.8092, 
  kr_sd = 180.5162, 
  kf_fumarase = 723.0053, 
  kr_fumarase = 395.3746,
  kf_md = 132.9511, 
  kr_md = 1557.9368, 
  kf_pta = 22.1455, 
  kr_pta = 237.7203,
  kf_ack = 4791.0944,
  kr_ack = 852.687,
  
  
  
  conserved_coa = 2.2090,
  conserved_quantity2 = 5.0533,
  
  k = 100,
  alpha = 0.0545/2.8
  
)


y0 <- c(  pyruvate = 3.66,
          oxaloacetate = 0.000487,
          acetyl_coa = 0.606,
          citrate = 1.96,
          cis_aconitate = 0.0161,
          isocitrate = 0.0367,
          aketo = 0.443,
          succinyl_coa = 0.233,
          succinate = 0.569,
          fumarate = 0.115,
          acetyl_p = 1.07,
          nad = 2.8,
          nadh = 0.0545)


ode_func <- function(t, y, parms) {
  
  with(as.list(c(y, parms)), {
    
    ## nadh = 0.0545
    atp = 8.13
    fadh2 = 0.0173
    co2 = 0.0752
    
    adp = 0.704
    fad = 0.173
    ## nad = 2.80
    pep = 0.231
    acetate = 0.005
    pi = 23.9
    
    pyruvate <- y[1]
    oxaloacetate <- y[2]
    acetyl_coa <- y[3]
    citrate <- y[4]
    cis_aconitate <- y[5]
    isocitrate <- y[6]
    aketo <- y[7]
    succinyl_coa <- y[8]
    succinate <- y[9]
    fumarate <- y[10]
    acetyl_p <- y[11]
    nad <- y[12]
    nadh <- y[13]
    
    
    coa_sh <- conserved_coa - acetyl_coa - succinyl_coa
    malate <- conserved_quantity2 - (oxaloacetate + citrate + cis_aconitate 
    + isocitrate + aketo + succinyl_coa + succinate + fumarate)
    
    
    
    # ---- reaction rates ----
    v_pdh <- compute_v_pdh(pyruvate, coa_sh, nad, acetyl_coa, nadh, co2,
                          km_co2_pdh, km_pyruvate_pdh, km_coa_pdh, km_nad_pdh, km_acetyl_coa_pdh, km_nadh_pdh,
                          enzyme_pdh, kf_pdh, kr_pdh, ki_nadh_pdh)
    
    
    v_cs <- compute_v_cs(oxaloacetate, acetyl_coa, citrate, coa_sh, nadh,
                        km_oxaloacetate_cs, km_acetyl_coa_cs, km_citrate_cs, km_coa_sh_cs,
                        enzyme_cs, kf_cs, kr_cs, ki_nadh_cs)
    
    
    v_aco <- compute_v_aco(citrate, cis_aconitate, isocitrate,
                                     km_citrate_aco, km_cis_aconitate_aco1, km_cis_aconitate_aco2, km_isocitrate_aco,
                                     enzyme_aco, kf1_aco, kr1_aco, kf2_aco, kr2_aco)
    
    v1_aco <- v_aco[1]
    v2_aco <- v_aco[2]
    
    
    v_id <- compute_v_id(isocitrate, nad, nadh, aketo, co2,
                        km_co2_id, km_isocitrate_id, km_nadh_id, km_nad_id, km_aketo_id,
                        enzyme_id, kf_id, kr_id)
    
    
    v_akd <- compute_v_akd(aketo, nad, coa_sh, succinyl_coa, nadh, co2,
                          km_co2_akd, km_aketo_akd, km_nad_akd, km_nadh_akd, km_coa_sh_akd, km_succinyl_coa_akd,
                          enzyme_akd, kf_akd, kr_akd)
    
    
    v_scs <- compute_v_scs(succinyl_coa, adp, succinate, coa_sh, atp, pi,
                          km_succinyl_coa_scs, km_adp_scs, km_succinate_scs, km_coa_sh_scs, km_atp_scs, km_pi_scs,
                          enzyme_scs, kf_scs, kr_scs)
    
    
    
    v_sd <- compute_v_sd(succinate, fad, fumarate, fadh2,
                        km_succinate_sd, km_fad_sd, km_fumarate_sd, km_fadh2_sd,
                        enzyme_sd, kf_sd, kr_sd)
    
    
    
    v_fumarase <- compute_v_fumarase(fumarate, malate,
                                    km_fumarate_fumarase, km_malate_fumarase, 
                                    enzyme_fumarase, kf_fumarase, kr_fumarase)
    
    
    
    v_md <- compute_v_md(malate, nad, oxaloacetate, nadh,
                        km_malate_md, km_nad_md, km_oxaloacetate_md, km_nadh_md,
                        enzyme_md, kf_md, kr_md)
    
    v_pk <- compute_v_pk(km_pyruvate_pk, km_atp_pk, km_pep_pk, km_adp_pk, kf_pk, kr_pk, enzyme_pk, pyruvate, atp, pep, adp, ki_atp_pk)
    
    
    v_ack <- compute_v_ack(acetyl_p, adp, acetate, atp, km_acetyl_p_ack, km_adp_ack, km_acetate_ack, km_atp_ack, kf_ack, kr_ack, enzyme_ack)
    
    
    v_pta <- compute_v_pta(acetyl_coa, pi, acetyl_p, coa_sh, km_acetyl_coa_pta, km_pi_pta, km_acetyl_p_pta, km_coa_pta, kf_pta, kr_pta, enzyme_pta)
    
    v_redox <- compute_v_redox(nadh, nad, k, alpha)
    
    # ODE system
    dPyr       <- v_pk - v_pdh
    dOAA       <- v_md - v_cs
    dAcCoA     <- v_pdh - v_cs - v_pta
    dCIT       <- v_cs - v1_aco
    dCIS       <- v1_aco - v2_aco
    dISO       <- v2_aco - v_id
    dAKG       <- v_id - v_akd
    dSucCoA    <- v_akd - v_scs
    dSuccinate <- v_scs - v_sd
    dFumarate  <- v_sd - v_fumarase
    dAcetylP   <- v_pta - v_ack
    dNAD       <- - v_pdh - v_id - v_akd - v_md - v_redox
    dNADH      <- v_pdh + v_id + v_akd + v_md + v_redox
    
    
    list(c(dPyr, dOAA, dAcCoA, dCIT, dCIS, dISO, dAKG, dSucCoA, dSuccinate, dFumarate, dAcetylP, dNAD, dNADH))
  })
  
  
  
  
}



t_obs <- seq(0, 60*60*2400, length.out = 2)

out <- ode(y = y0, times = t_obs, func = ode_func, parms = params_true,
           method = "vode", rtol = 1e-6, atol = 1e-9)



y_true <- out[nrow(out), -1]

y_true


v_ack <- compute_v_ack(y_true[11], 0.704, 0.005, 8.13, params_true$km_acetyl_p_ack, params_true$km_adp_ack, params_true$km_acetate_ack, params_true$km_atp_ack, params_true$kf_ack, params_true$kr_ack, params_true$enzyme_ack)


v_ack


v_sd <- compute_v_sd(y_true[9], 0.173, y_true[10], 0.0173,
                     params_true$km_succinate_sd, params_true$km_fad_sd, params_true$km_fumarate_sd, params_true$km_fadh2_sd,
                     params_true$enzyme_sd, params_true$kf_sd, params_true$kr_sd)

v_sd


set.seed(123)
sigma <- 0.15
log10_y_obs <- rnorm(length(y_true), log10(y_true), sigma)

y_obs <- 10.^log10_y_obs

y_obs

log10_v_acka_obs <-  rnorm(1, log10(v_ack), sigma)
  
v_acka_obs <- 10.^log10_v_acka_obs

v_acka_obs



log10_v_sd_obs <-  rnorm(1, log10(v_sd), sigma)

v_sd_obs <- 10.^log10_v_sd_obs

v_sd_obs



