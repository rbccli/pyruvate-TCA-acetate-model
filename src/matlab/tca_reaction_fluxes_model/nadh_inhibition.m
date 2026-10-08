function inh_nadh = nadh_inhibition(nadh_con, ki_nadh)
    inh_nadh = 1 ./ (1 + nadh_con ./ ki_nadh);
    %inh_nadh = 1;
end