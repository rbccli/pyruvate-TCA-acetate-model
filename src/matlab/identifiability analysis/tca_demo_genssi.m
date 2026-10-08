function model = tca_demo_genssi()


    syms pyruvate acetyl_coa oxaloacetate 
    syms citrate cis_aconitrate isocitrate alpha_ketoglutarate succinyl_coa succinate fumarate acetyl_p nad nadh

    syms enzyme_pdh enzyme_pk enzyme_cs
    syms enzyme_aco enzyme_akd enzyme_id
    syms enzyme_scs enzyme_sd enzyme_pta 
    syms enzyme_fumarase enzyme_md k alp enzyme_ack %qh2 q acetate
    
    %qh2 q enzyme_ack 
    
    % k alp
    


    model.sym.x = [pyruvate; acetyl_coa; oxaloacetate; 
        citrate; cis_aconitrate; isocitrate; alpha_ketoglutarate; succinyl_coa; succinate; fumarate; acetyl_p; nad; nadh];

    model.sym.p = [ enzyme_pdh; enzyme_pk; enzyme_aco; enzyme_cs; enzyme_id;
        enzyme_akd; enzyme_scs; enzyme_sd; enzyme_fumarase; enzyme_md; enzyme_pta; k; alp; enzyme_ack]; 
        % qh2; q; enzyme_ack; k; alp; acetate];
    

    conserved_coa = 2.2;
    conserved_quantity2 = 5.05;

    %enzyme_ack = 0.1;

    % enzyme_fumarase = 0.1;

    % k = 100000;
    % 
    % alp = 0.0545/2.8;

    %% parameter

    acetate = 1;

    % nadh = 0.0545; %% lower bound value
    atp = 8.13; %% fake data, true: 9.63 lb: 8.13 
    qh2 = 0.0173; %% qh2 no data 
    co2 = 0.0752; 

    adp = 0.704;
    q = 0.173; %% q no data
    % nad = 2.80;
    pep = 0.231;
    pi = 23.9;


    ki_nadh_pdh = 0.001;
    ki_nadh_cs = 0.0028;
    ki_atp_pk = 4;


    % michaelis constant
    % pk
    km_pyruvate_pk = 0.0505;
    km_pep_pk = 0.2666;
    km_atp_pk = 16.1752;
    km_adp_pk = 0.1919;


    % pdh
    km_pyruvate_pdh = 0.3296; 
    km_coa_pdh = 0.0471; 
    km_nad_pdh = 0.0619; 
    km_nadh_pdh = 0.215;
    km_acetyl_coa_pdh = 0.215;
    km_co2_pdh = 0.215;

    % pta
    km_acetyl_coa_pta = 0.0491;
    km_pi_pta = 1.9606;
    km_acetyl_p_pta = 0.6297;
    km_coa_pta = 0.0521;
    
    % ack
    km_acetyl_p_ack = 0.0938;
    km_adp_ack = 0.2799;
    km_acetate_ack = 7.9678;
    km_atp_ack = 0.096;
    

    % cs
    km_citrate_cs = 0.0018; 
    km_coa_cs = 44914.6131; % 0.0027; 
    km_oxaloacetate_cs = 0.0119; 
    km_acetyl_coa_cs = 0.0521; 


    % aconitase
    km_citrate_aco = 0.6464; 
    km_cis_aconitate_aco1 = 0.0095;
    km_cis_aconitate_aco2 = 0.2098;
    km_isocitrate_aco = 0.0872; 

    % isocitrate dehydrogenase
    km_isocitrate_id = 0.0109; 
    km_nad_id = 1.5211; 
    km_aketo_id = 0.4744; 
    km_nadh_id = 0.007; % 24.3; 
    km_co2_id = 2.4317;


    % akd checked not much different
    km_nad_akd = 0.0512; 
    km_a_ketoglutarate_akd = 0.0026; 
    km_coa_akd = 0.0512; 
    km_nadh_akd = 0.1953; 
    km_succinyl_coa_akd = 0.1953; 
    km_co2_akd = 0.1953;


    % scs finished check 
    km_coa_sh_scs = 2.0659; 
    km_succinyl_coa_scs = 0.0316; 
    km_adp_scs = 0.3244; 
    km_succinate_scs = 0.1444; 
    km_atp_scs = 0.0426; 
    km_pi_scs = 0.8953;


    % sd finished check 
    km_succinate_sd = 0.0021;  
    km_q_sd = 0.002; % q
    km_fumarate_sd = 0.0233; 
    km_qh2_sd = 12.5427; % qh2


    % fumarase
    km_fumarate_famarase = 0.4222; 
    km_malate_fumarase = 0.6905; 

    % malate dehydrogenase dont change same ecoli exp
    km_nad_md = 0.418; 
    km_oxaloacetate_md = 0.0244; 
    km_malate_md = 3.8086; 
    km_nadh_md = 0.0301; 





    % forward and reverse k_cat 

    kf_pk = 555.4713;
    kr_pk = 0.4124;
    kf_pdh = 36.1192; 
    kr_pdh = 0;
    kf_cs = 141.5176; 
    kr_cs = 5.8632; 
    kf1_aco = 5.7836; 
    kf2_aco = 4.5447; 
    kr1_aco = 3.8424; 
    kr2_aco = 1.0233; 
    kf_id = 82.0188; 
    kr_id = 14.7374;
    kf_akd = 2.0345; 
    kr_akd = 0.0011;
    kf_scs = 35.8092; 
    kr_scs = 60.5031; 
    kf_sd = 32.8092; 
    kr_sd = 180.5162; 
    kf_fumarase = 723.0053; 
    kr_fumarase = 395.3746;
    kf_md = 132.9511; 
    kr_md = 1557.9368; 
    kf_pta = 22.1455; 
    kr_pta = 237.7203;
    kf_ack = 4791.0944;
    kr_ack = 852.687;

    
    
    coa = conserved_coa - acetyl_coa - succinyl_coa;
    malate = conserved_quantity2 - oxaloacetate - citrate - cis_aconitrate - isocitrate - alpha_ketoglutarate - succinyl_coa - succinate - fumarate;

    
    %% enzyme pdh
    % pyruvate + nad + coa -> acetyl_coa + nadh + co2

    pyruvate_tilde = pyruvate / km_pyruvate_pdh;
    coa_tilde = coa / km_coa_pdh;
    nad_tilde = nad / km_nad_pdh;
    acetyl_coa_tilde = acetyl_coa/km_acetyl_coa_pdh;
    nadh_tilde = nadh / km_nadh_pdh;
    co2_tilde = co2/ km_co2_pdh;
    
    denominator1 = (1 + pyruvate_tilde) .* (1 + coa_tilde) .* (1 + nad_tilde) + (1 + acetyl_coa_tilde) .* (1 + nadh_tilde) .* (1 + co2_tilde) - 1;
    
    v1 = (kf_pdh .* pyruvate_tilde .* coa_tilde .* nad_tilde - kr_pdh .* acetyl_coa_tilde .* nadh_tilde .* co2_tilde)./ denominator1 * enzyme_pdh .* 1 ./ (1 + nadh ./ ki_nadh_pdh);


    %% enzyme pk
    % pep + adp -> pyruvate + atp

    pyruvate_tilde = pyruvate ./ km_pyruvate_pk;
    pep_tilde = pep ./ km_pep_pk;
    atp_tilde = atp ./ km_atp_pk;
    adp_tilde = adp ./ km_adp_pk;

    denominator2 = (1+pyruvate_tilde).*(1+atp_tilde)+(1+pep_tilde).*(1+adp_tilde)-1;

    v2 = (kf_pk .* pep_tilde .* adp_tilde - kr_pk .* pyruvate_tilde .* atp_tilde) ./ denominator2 * enzyme_pk .* 1 ./ (1 + atp ./ ki_atp_pk);


    %% enzyme cs
    % oxalocetate + acetyle-coa + h20 <-> citrate + coa-sh

    ox_tilde = oxaloacetate / km_oxaloacetate_cs;
    acc_tilde = acetyl_coa / km_acetyl_coa_cs;
    cit_tilde = citrate / km_citrate_cs;
    coa_tilde = coa / km_coa_cs;

    denominator3 = (1 + ox_tilde) .* (1 + acc_tilde) + (1 + cit_tilde) .* (1 + coa_tilde) - 1;

    v3 = (kf_cs * ox_tilde .* acc_tilde - kr_cs * cit_tilde .* coa_tilde) ./ denominator3 * enzyme_cs .* 1 ./ (1 + nadh ./ ki_nadh_cs);



    %% enzyme aco
    % citrate <-> cis-aconitrate + h20 <-> isocitrate

    cit_tilde = citrate / km_citrate_aco;
    cis_tilde1 = cis_aconitrate / km_cis_aconitate_aco1;
    cis_tilde2 = cis_aconitrate / km_cis_aconitate_aco2;
    iso_tilde = isocitrate / km_isocitrate_aco;

    v4 = (kf1_aco * cit_tilde - kr1_aco * cis_tilde1) ./ ((1 + cit_tilde) + (1 + cis_tilde1) - 1) * enzyme_aco;
    v5 = (kf2_aco * cis_tilde2 - kr2_aco * iso_tilde) ./ ((1 + cis_tilde2) + (1 + iso_tilde) - 1) * enzyme_aco;



    %% enzyme id
    % isocitrate + nad -> nadh + alpha-ketoglutarate + co2 + h+

    iso_tilde = isocitrate / km_isocitrate_id;
    nad_tilde = nad / km_nad_id;
    nadh_tilde = nadh / km_nadh_id;
    aketo_tilde = alpha_ketoglutarate/km_aketo_id;
    co2_tilde = co2/ km_co2_id;


    denominator6 = (1 + iso_tilde) .* (1 + nad_tilde) + (1 + nadh_tilde) .* (1 + aketo_tilde) .* (1 + co2_tilde)- 1;

    v6 = (kf_id .* iso_tilde .* nad_tilde - kr_id .* nadh_tilde .* aketo_tilde .* co2_tilde)./ denominator6 * enzyme_id;


    
    %% enzyme akd
    % alpha-ketoglutarate + nad + coa-sh -> succinyl-coa + nadh + co2 + h+

    ak_tilde = alpha_ketoglutarate / km_a_ketoglutarate_akd;
    nad_tilde = nad / km_nad_akd;
    coa_tilde = coa / km_coa_akd;
    succinyl_coa_tilde = succinyl_coa / km_succinyl_coa_akd;
    nadh_tilde = nadh / km_nadh_akd;
    co2_tilde = co2/ km_co2_akd;

    denominator7 = (1 + nad_tilde) .* (1 + ak_tilde) .* (1 + coa_tilde) + (1 + succinyl_coa_tilde) .* (1 + nadh_tilde) .* (1 + co2_tilde) - 1;

    v7 = (kf_akd * ak_tilde .* coa_tilde .* nad_tilde - kr_akd * succinyl_coa_tilde .* nadh_tilde .* co2_tilde)./ denominator7 * enzyme_akd;

    %% enzyme scs
    % succinyl-coa + adp <-> succinate + coa + atp

    sucCoA_tilde = succinyl_coa / km_succinyl_coa_scs;
    adp_tilde = adp / km_adp_scs;
    suc_tilde = succinate / km_succinate_scs;
    coa_tilde = coa / km_coa_sh_scs;
    atp_tilde = atp / km_atp_scs;
    pi_tilde = pi / km_pi_scs;
    
    denominator8 = (1 + sucCoA_tilde) .* (1 + adp_tilde) + (1 + suc_tilde) .* (1 + coa_tilde) .* (1 + atp_tilde) .* (1 + pi_tilde) - 1;
    
    v8 = (kf_scs * sucCoA_tilde .* adp_tilde - kr_scs * suc_tilde .* coa_tilde .* atp_tilde .* pi_tilde) ./ denominator8 * enzyme_scs;



    %% enzyme sd
    % succinate + fad <-> fumarate + fadh2

    suc_tilde = succinate / km_succinate_sd;
    uq_tilde = q / km_q_sd;
    fum_tilde = fumarate / km_fumarate_sd;
    uqh2_tilde = qh2 / km_qh2_sd;
    
    denominator9 = (1 + suc_tilde) .* (1 + uq_tilde) + (1 + fum_tilde) .* (1 + uqh2_tilde) - 1;

    v9 = (kf_sd * suc_tilde .* uq_tilde - kr_sd * fum_tilde .* uqh2_tilde) ./ denominator9 * enzyme_sd;


    
    %% enzyme fumarase
    % fumarate + h20 <-> malate

    fum_tilde = fumarate / km_fumarate_famarase;
    mal_tilde = malate / km_malate_fumarase;

    denominator10 = (1 + fum_tilde) + (1 + mal_tilde) - 1;

    v10 = (kf_fumarase * fum_tilde - kr_fumarase * mal_tilde) ./ denominator10 * enzyme_fumarase;



    %% enzyme malate dehydrogenase

    % malate + nad <-> oxaloacetate + nadh

    mal_tilde = malate / km_malate_md;
    nad_tilde = nad / km_nad_md;
    oaa_tilde = oxaloacetate / km_oxaloacetate_md;
    nadh_tilde = nadh / km_nadh_md;
    
    denominator11 = (1 + mal_tilde) .* (1 + nad_tilde) + (1 + oaa_tilde) .* (1 + nadh_tilde) - 1;

    v11 = (kf_md * mal_tilde .* nad_tilde - kr_md * oaa_tilde .* nadh_tilde) ./ denominator11 * enzyme_md;


   
    %% enzyme pta

    acetyl_coa_tilde = acetyl_coa / km_acetyl_coa_pta;
    pi_tilde = pi / km_pi_pta;
    acetyl_p_tilde = acetyl_p / km_acetyl_p_pta;
    coa_tilde = coa / km_coa_pta;


    denominator12 = (1 + acetyl_coa_tilde).* (1 + pi_tilde) + (1 + acetyl_p_tilde) .* (1 + coa_tilde) - 1;

    v12 = (kf_pta .* acetyl_coa_tilde .* pi_tilde - kr_pta .* acetyl_p_tilde .* coa_tilde)./ denominator12 * enzyme_pta;




    %% enzyme ack
    acetyl_p_tilde = acetyl_p / km_acetyl_p_ack;
    adp_tilde = adp / km_adp_ack;

    atp_tilde = atp / km_atp_ack;
    acetate_tilde = acetate / km_acetate_ack;


    denominator13 = (1 + acetyl_p_tilde) .* (1 + adp_tilde) + (1 + acetate_tilde) .* (1 + atp_tilde) - 1;

    v13 = (kf_ack .* acetyl_p_tilde .* adp_tilde - kr_ack .* acetate_tilde .* atp_tilde)./ denominator13 * enzyme_ack;



    %% redox reaction

    v_redox = k * alp * nad - k * nadh;
    
   

    %% model

    % model.sym.g = sym(zeros(13, 1));

    model.sym.g = [];

    model.sym.xdot = [

        -v1 + v2 % pyruvate
        v1 - v3 - v12 % acetyl-coa
        -v3 + v11 % oxaloacetate
        v3 - v4 % citrate
        v4 - v5 % cis_aconitate
        v5 - v6 % isocitrate
        v6 - v7 % alpha_ketoglutarate
        v7 - v8 % succinyl_coa
        v8 - v9 % succinate
        v9 - v10 % fumarase
        v12 - v13 % acetyl_p
        - v1 - v6 - v7 - v11 - v_redox % nad+
        v1 + v6 + v7 + v11 + v_redox % nadh%
        ];

    model.sym.x0 = [3.66; 0.000487; 0.606; 1.96; 0.0161; 0.0367; 0.443; 0.233; 0.569; 0.115; 1.07; 2.8; 0.0545];

    model.sym.y = [pyruvate; acetyl_coa; oxaloacetate; citrate; cis_aconitrate; isocitrate; alpha_ketoglutarate; succinyl_coa; succinate; fumarate; acetyl_p; nadh; nad; v3; v4; v5; v6; v7; v8; v9; v10; v11; v12; v13];


end