function v = compute_v_scs(succinyl_coa_con, adp_con, succinate_con, coa_sh_con, atp_con, pi_con, ...
km_succinyl, km_adp, km_succinate, km_coa, km_atp, km_pi, enzyme_scs, kf_scs, kr_scs)
% Reaction 6: succinyl-coa + adp + pi <-> succinate + coa-sh + atp

    sucCoA_tilde = succinyl_coa_con / km_succinyl;
    adp_tilde = adp_con / km_adp;
    suc_tilde = succinate_con / km_succinate;
    coa_tilde = coa_sh_con / km_coa;
    atp_tilde = atp_con / km_atp;
    pi_tilde = pi_con / km_pi;

    
    denominator = (1 + sucCoA_tilde) .* (1 + adp_tilde) .* (1 + pi_tilde)+ (1 + suc_tilde) .* (1 + coa_tilde) .* (1 + atp_tilde) - 1;
    
    v = (kf_scs * sucCoA_tilde .* adp_tilde .* pi_tilde - kr_scs * suc_tilde .* coa_tilde .* atp_tilde) ./ denominator * enzyme_scs;


end