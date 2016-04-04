function roots_vec = find_roots_binsearch_fzero(func, left_bound, right_bound, mid)

mid = fzero(func, mid);

if ( ~isnan(mid) )
    left_root_vec = find_roots_binsearch_fzero(func, left_bound, mid-eps(mid));
    right_root_vec = find_roots_binsearch_fzero(func, mid+eps(mid));
    roots_vec = [left_root_vec, mid, right_root_vec];
else
    roots_vec = [];
end



if ( ~keep_going )
    roots_vec = [ known_roots_vec, left_root];
else
    mid = 0.5*( left_bound + right_bound);
    %mid_root = vpasolve(vpa_func == 0, vpa_var, mid)
    
    % [FIX] Change to key off whether fzero found a root in the interval or
    % not. Use bounds for solutions rather than initial point
    mid_root = fzero(func, mid);
    fprintf('%.3g  %.3g  %.3g  ::  %.3g  %.3g  %.3g  ::  %d  %.3g\n', left_bound, mid, right_bound, left_root, mid_root, right_root, keep_going, raw_err);
    roots_vec = [known_roots_vec, mid_root];
    lower_roots_vec = find_roots_binsearch_fzero(func, left_bound, mid, left_root, mid_root, raw_err_tol, roots_vec);
    upper_roots_vec = find_roots_binsearch_fzero(func, mid, right_bound, mid_root, right_root, raw_err_tol, roots_vec);
    
    roots_vec = [lower_roots_vec, roots_vec, upper_roots_vec];
end