function inh_atp = atp_inhibition(atp_con, ki_atp)
    inh_atp = 1 ./ (1 + atp_con ./ ki_atp);
    inh_atp = 1;
end
