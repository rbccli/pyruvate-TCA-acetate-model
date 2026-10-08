function v = compute_v_akd(a_keto_con, nad_con, coa_sh_con, succinyl_coa_con, nadh_con, co2_con, km_co2, km_ak, km_nad, km_nadh, km_coa, km_succinyl_coa, enzyme_akd, kf_akd, kr_akd)
% Reaction 5: alpha-ketoglutarate + nad + coa-sh -> succinyl-coa + nadh + co2 + h+

    ak_tilde = a_keto_con / km_ak;
    nad_tilde = nad_con / km_nad;
    coa_tilde = coa_sh_con / km_coa;
    succinyl_coa_tilde = succinyl_coa_con / km_succinyl_coa;
    nadh_tilde = nadh_con / km_nadh;
    co2_tilde = co2_con/ km_co2;

    denominator = (1 + nad_tilde) .* (1 + ak_tilde) .* (1 + coa_tilde) + (1 + succinyl_coa_tilde) .* (1 + nadh_tilde) .* (1 + co2_tilde) - 1;

    v = (kf_akd * ak_tilde .* coa_tilde .* nad_tilde - kr_akd .* succinyl_coa_tilde .* nadh_tilde .* co2_tilde)./ denominator * enzyme_akd;

end
