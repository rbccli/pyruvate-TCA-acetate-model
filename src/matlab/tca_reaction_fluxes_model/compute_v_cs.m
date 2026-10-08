function v = compute_v_cs(oxaloacetate_con, acetyle_coa_con, citrate_con, coa_sh_con, nadh_con, ...
km_ox, km_acc, km_cit, km_coa_sh, enzyme_cs, kf_cs, kr_cs, ki_nadh, baseline)
% Reaction 2: oxalocetate + acetyle-coa + h20 <-> citrate + coa-sh


    ox_tilde = oxaloacetate_con / km_ox;
    acc_tilde = acetyle_coa_con / km_acc;
    cit_tilde = citrate_con / km_cit;
    coa_tilde = coa_sh_con / km_coa_sh;

    denominator = (1 + ox_tilde) .* (1 + acc_tilde) + (1 + cit_tilde) .* (1 + coa_tilde) - 1;

    v = (kf_cs * ox_tilde .* acc_tilde - kr_cs * cit_tilde .* coa_tilde) ./ denominator * enzyme_cs .* nadh_inhibition(nadh_con, ki_nadh);


    % 0.0488656195462478; %
    % baseline; %
end
