function v = compute_v_fumarase(fumarate_con, malate_con, km_fum, km_mal, enzyme_fum, kf_fum, kr_fum)
% Reaction 8: fumarate + h20 <-> malate

    fum_tilde = fumarate_con / km_fum;
    mal_tilde = malate_con / km_mal;

    denominator = (1 + fum_tilde) + (1 + mal_tilde) - 1;

    v = (kf_fum * fum_tilde - kr_fum * mal_tilde) ./ denominator * enzyme_fum;

end
