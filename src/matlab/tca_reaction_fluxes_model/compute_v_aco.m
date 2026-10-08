function [v1, v2] = compute_v_aco(citrate_con, cis_aconitate_con, isocitrate_con, ...
km_cit, km_cis1, km_cis2, km_iso, enzyme_aco, kf1, kr1, kf2, kr2)
% Reaction 3: citrate <-> cis-aconitrate + h20 <-> isocitrate

    cit_tilde = citrate_con / km_cit;
    cis_tilde1 = cis_aconitate_con / km_cis1;
    
    cis_tilde2 = cis_aconitate_con / km_cis2;
    iso_tilde = isocitrate_con / km_iso;
    
    v1 = (kf1 * cit_tilde - kr1 * cis_tilde1) ./ ((1 + cit_tilde) + (1 + cis_tilde1) - 1) * enzyme_aco;
    v2 = (kf2 * cis_tilde2 - kr2 * iso_tilde) ./ ((1 + cis_tilde2) + (1 + iso_tilde) - 1) * enzyme_aco;

    
end