function v = compute_v_md(malate_con, nad_con, oaa_con, nadh_con, ...
km_mal, km_nad, km_oaa, km_nadh, enzyme_md, kf_md, kr_md)
% Reaction 9: malate + nad <-> oxaloacetate + nadh + h+

    mal_tilde = malate_con / km_mal;
    nad_tilde = nad_con / km_nad;
    oaa_tilde = oaa_con / km_oaa;
    nadh_tilde = nadh_con / km_nadh;
    
    denominator = (1 + mal_tilde) .* (1 + nad_tilde) + (1 + oaa_tilde) .* (1 + nadh_tilde) - 1;

    v = (kf_md * mal_tilde .* nad_tilde - kr_md * oaa_tilde .* nadh_tilde) ./ denominator * enzyme_md;


end