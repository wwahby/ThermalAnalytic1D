function roots_vec = find_roots_binsearch(vpa_func, vpa_var, left_bound, right_bound, left_root, right_root, raw_err_tol, known_roots_vec)


%root_left = vpasolve(vpa_func == 0, vpa_var, min_bound);
%root_right = vpasolve(vpa_func == 0, vpa_var, max_bound);
raw_err = abs(double(left_root) - double(right_root));
keep_going = (raw_err > raw_err_tol) % keep going if the roots are different (to whatever tolerance)
if ( ~keep_going )
    roots_vec = [ known_roots_vec, left_root];
    return;
else
    mid = 1/2* ( left_bound + right_bound)
    mid_root = vpasolve(vpa_func == 0, vpa_var, [mid inf])
    roots_vec = [known_roots_vec, mid_root]
    lower_roots_vec = find_roots_binsearch(vpa_func, vpa_var, left_bound, mid, left_root, mid_root, raw_err_tol, roots_vec);
    upper_roots_vec = find_roots_binsearch(vpa_func, vpa_var, mid, right_bound, mid_root, right_root, raw_err_tol, roots_vec);
    
    roots_vec = [lower_roots_vec, roots_vec, upper_roots_vec];
end