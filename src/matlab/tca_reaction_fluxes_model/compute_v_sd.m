function v = compute_v_sd(succinate_con, uq_con, fumarate_con, uqh2_con, ...
km_suc, km_uq, km_fum, km_uqh2, enzyme_sd, kf_sd, kr_sd)
% Reaction 7: succinate + fad <-> fumarate + fadh2

    suc_tilde = succinate_con / km_suc;
    uq_tilde = uq_con / km_uq;
    fum_tilde = fumarate_con / km_fum;
    uqh2_tilde = uqh2_con / km_uqh2;
    
    denominator = (1 + suc_tilde) .* (1 + uq_tilde) + (1 + fum_tilde) .* (1 + uqh2_tilde) - 1;

    v = (kf_sd * suc_tilde .* uq_tilde - kr_sd * fum_tilde .* uqh2_tilde) ./ denominator * enzyme_sd;

    
end