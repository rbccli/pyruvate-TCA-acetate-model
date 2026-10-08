function p_full = insert_params(p_fixed, p_free, idx_free)
    p_full = p_fixed;
    p_full(idx_free) = p_free;
end