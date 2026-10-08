function [dydt] = function_fsolve(y, params)


    % mitobolite defining
    pyruvate = y(1);
    oxaloacetate = y(2);
    acetyl_coa = y(3);
    citrate = y(4);
    cis_aconitate = y(5);
    isocitrate = y(6);
    aketo = y(7);
    succinyl_coa = y(8);
    succinate = y(9);
    fumarate = y(10);
    acetyl_p = y(11);
   
    coa_sh = params.conserved_coa - acetyl_coa - succinyl_coa;
    malate = params.conserved_quantity2 - oxaloacetate - citrate - cis_aconitate - isocitrate - aketo - succinyl_coa - succinate - fumarate;

    % malate = y(14); %%SANITY CHECK ONLY

    
    nad = y(12); % params.nad_pool * (1 - params.nadh_ratio);
    nadh = y(13); % params.nad_pool * params.nadh_ratio;


    acetate = params.acetate_constant;
    pi = params.pi_constant;
    atp = params.adp_pool * params.atp_ratio;
    qh2 = params.q_pool * params.qh2_ratio;
    adp = params.adp_pool * (1 - params.atp_ratio);
    q = params.q_pool * (1 - params.qh2_ratio);


    co2 = params.co2_constant;
    pep = params.pep_constant;

    
    % reaction rates
    % v_pdh = compute_v_pdh(pyruvate, coa_sh, nad, acetyl_coa, nadh, co2, ...
    %     params.km_co2_pdh, params.km_pyruvate_pdh, params.km_coa_pdh, params.km_nad_pdh, params.km_acetyl_coa_pdh, params.km_nadh_pdh, ...
    %     params.enzyme_pdh, params.kf_pdh, params.kr_pdh, params.ki_nadh_pdh);

    v_pdh = compute_v_pdh_beta(pyruvate, coa_sh, nad, acetyl_coa, nadh, co2, ...
    params.km_co2_pdh, params.km_pyruvate_pdh, params.km_coa_pdh, params.km_nad_pdh, params.km_acetyl_coa_pdh, params.km_nadh_pdh, ...
    params.enzyme_pdh, params.kf_pdh, params.kr_pdh, params.ki_nadh_pdh, params.beta, params.baseline_pdh);

    v_cs = compute_v_cs(oxaloacetate, acetyl_coa, citrate, coa_sh, nadh, ...
        params.km_oxalocetate_cs, params.km_acetyl_coa_cs, params.km_citrate_cs, params.km_coa_sh_cs, ...
        params.enzyme_cs, params.kf_cs, params.kr_cs, params.ki_nadh_cs, params.baseline_cs);
    

    [v1_aco, v2_aco] = compute_v_aco(citrate, cis_aconitate, isocitrate, ...
        params.km_citrate_aco, params.km_cis_aconitate_aco1, params.km_cis_aconitate_aco2, params.km_isocitrate_aco, ...
        params.enzyme_aco, params.kf1_aco, params.kr1_aco, params.kf2_aco, params.kr2_aco);

    
    v_id = compute_v_id(isocitrate, nad, nadh, aketo, co2, ...
        params.km_co2_id, params.km_isocitrate_id, params.km_nadh_id, params.km_nad_id, params.km_aketo_id, ...
        params.enzyme_id, params.kf_id, params.kr_id);


    v_akd = compute_v_akd(aketo, nad, coa_sh, succinyl_coa, nadh, co2, ...
        params.km_co2_akd, params.km_aketo_akd, params.km_nad_akd, params.km_nadh_akd, params.km_coa_sh_akd, params.km_succinyl_coa_akd, ...
        params.enzyme_akd, params.kf_akd, params.kr_akd);
    

    v_scs = compute_v_scs(succinyl_coa, adp, succinate, coa_sh, atp, pi, ...
        params.km_succinyl_coa_scs, params.km_adp_scs, params.km_succinate_scs, params.km_coa_sh_scs, params.km_atp_scs, params.km_pi_scs,...
        params.enzyme_scs, params.kf_scs, params.kr_scs);

    
    v_sd = compute_v_sd(succinate, q, fumarate, qh2, ...
        params.km_succinate_sd, params.km_q_sd, params.km_fumarate_sd, params.km_qh2_sd, ...
        params.enzyme_sd, params.kf_sd, params.kr_sd);
    
 
    v_fumarase = compute_v_fumarase(fumarate, malate, ...
        params.km_fumarate_famarase, params.km_malate_fumarase, ...
        params.enzyme_fumarase, params.kf_fumarase, params.kr_fumarase);

    
    v_md = compute_v_md(malate, nad, oxaloacetate, nadh, ...
        params.km_malate_md, params.km_nad_md, params.km_oxalocetate_md, params.km_nadh_md, ...
        params.enzyme_md, params.kf_md, params.kr_md);

   
    v_pk = compute_v_pk(params.km_pyruvate_pk, params.km_atp_pk, params.km_pep_pk, params.km_adp_pk, params.kf_pk, params.kr_pk, params.enzyme_pk, pyruvate, atp, pep, adp, params.ki_atp_pk);

    
    v_ack = compute_v_ack(acetyl_p, adp, acetate, atp, params.km_acetyl_p_ack, params.km_adp_ack, params.km_acetate_ack, params.km_atp_ack, params.kf_ack, params.kr_ack, params.enzyme_ack);


    v_pta = compute_v_pta(acetyl_coa, pi, acetyl_p, coa_sh, params.km_acetyl_coa_pta, params.km_pi_pta, params.km_acetyl_p_pta, params.km_coa_pta, params.kf_pta, params.kr_pta, params.enzyme_pta);


    v_redox = compute_v_redox(nadh, nad, params.k, params.alpha);


    dydt = zeros(13,1);
    
    dydt(1) = v_pk - v_pdh; % pyruvate
    dydt(2) = v_md - v_cs; % oxalocetate
    dydt(3) = v_pdh - v_cs - v_pta; % act_coa
    dydt(4) = v_cs - v1_aco; % citrate
    dydt(5) = v1_aco - v2_aco; % cis_aconitate
    dydt(6) = v2_aco - v_id; % isocitrate
    dydt(7) = v_id - v_akd; % aketo
    dydt(8) = v_akd - v_scs; % succ_coa
    dydt(9) = v_scs - v_sd; % succ
    dydt(10) = v_sd - v_fumarase; % fumarate
    dydt(11) = v_pta - v_ack; % acetyl-p
    % dydt(14) = v_fumarase - v_md; % SANITY CHECK ONLY: MALATE
    dydt(12) = - v_pdh - v_id - v_akd - v_md - v_redox; % nad+
    dydt(13) = v_pdh + v_id + v_akd + v_md + v_redox; % nadh


end

