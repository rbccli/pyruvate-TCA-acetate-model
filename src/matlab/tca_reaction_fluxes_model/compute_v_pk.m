function v = compute_v_pk(km_pyruvate, km_atp, km_pep, km_adp, kf_pk, kr_pk, enzyme_pk, pyruvate_con, atp_con, pep_con, adp_con, ki_atp)

    % reaction: pyruvate + ATP → oxaloacetate + ADP
    
    pyruvate_tilde = pyruvate_con ./ km_pyruvate;
    pep_tilde = pep_con ./ km_pep;
    atp_tilde = atp_con ./ km_atp;
    adp_tilde = adp_con ./ km_adp;

    denominator = (1+pyruvate_tilde).*(1+atp_tilde)+(1+pep_tilde).*(1+adp_tilde)-1;


    v = (kf_pk .* pep_tilde .* adp_tilde - kr_pk .* atp_tilde .* pyruvate_tilde) ./ denominator * enzyme_pk .* atp_inhibition(atp_con, ki_atp);



    % 0.329760923330585; %
end