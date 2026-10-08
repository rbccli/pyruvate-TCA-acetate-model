function v = compute_v_id(isocitrate_con, nad_con, nadh_con, aketo_con, co2_con, km_co2, km_iso, km_nadh, km_nad, km_aketo, enzyme_id, kf_id, kr_id)
% Reaction 4: isocitrate + nad -> nadh + alpha-ketoglutarate + co2 + h+

    iso_tilde = isocitrate_con / km_iso;
    nad_tilde = nad_con / km_nad;
    nadh_tilde = nadh_con / km_nadh;
    aketo_tilde = aketo_con/km_aketo;
    co2_tilde = co2_con/ km_co2;


    denominator = (1 + iso_tilde) .* (1 + nad_tilde) + (1 + nadh_tilde) .* (1 + aketo_tilde) .* (1 + co2_tilde)- 1;

    v = (kf_id * iso_tilde .* nad_tilde - kr_id .* nadh_tilde .* aketo_tilde .* co2_tilde)./ denominator * enzyme_id;

end
