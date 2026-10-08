function v = compute_v_pta(acetyl_coa_con, pi_con, acetyl_p_con, coa_sh_con, km_acetyl_coa_pta, km_pi_pta, km_acetyl_p_pta, km_coa_pta, kf_pta, kr_pta, enzyme_pta)
% Reaction 5: alpha-ketoglutarate + nad + coa-sh -> succinyl-coa + nadh + co2

    acetyl_coa_tilde = acetyl_coa_con / km_acetyl_coa_pta;
    pi_tilde = pi_con / km_pi_pta;
    acetyl_p_tilde = acetyl_p_con / km_acetyl_p_pta;
    coa_tilde = coa_sh_con / km_coa_pta;


    denominator = (1 + acetyl_coa_tilde).* (1 + pi_tilde) + (1 + acetyl_p_tilde) .* (1 + coa_tilde) - 1;

    v = (kf_pta .* acetyl_coa_tilde .* pi_tilde - kr_pta .* acetyl_p_tilde .* coa_tilde)./ denominator * enzyme_pta;

end

