function v = compute_v_ack(acetyl_p_con, adp_con, acetate_con, atp_con, km_acetyl_p_ack, km_adp_ack, km_acetate_ack, km_atp_ack, kf_ack, kr_ack, enzyme_ack)
% Reaction 5: alpha-ketoglutarate + nad + coa-sh <-> succinyl-coa + nadh + co2 + h+


    acetyl_p_tilde = acetyl_p_con / km_acetyl_p_ack;
    adp_tilde = adp_con / km_adp_ack;

    atp_tilde = atp_con / km_atp_ack;
    acetate_tilde = acetate_con / km_acetate_ack;


    denominator = (1 + acetyl_p_tilde) .* (1 + adp_tilde) + (1 + acetate_tilde) .* (1 + atp_tilde) - 1;

    v = (kf_ack .* acetyl_p_tilde .* adp_tilde - kr_ack .* acetate_tilde .* atp_tilde)./ denominator * enzyme_ack;

end
