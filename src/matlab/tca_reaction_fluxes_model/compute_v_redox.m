function v = compute_v_redox(nadh, nad, k, alpha)

    % k+ = k;
    % k- = alpha * k;
    
    v = k * alpha* nad - k * nadh;

end
