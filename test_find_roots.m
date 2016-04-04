% test find_roots_binsearch

% syms x;
% 
% testfunc = vpa(sin(x));
% vpa_func = testfunc;
% vpa_var = x;
vpa_var = 1;

testfunc = @(x) sin(x);
vpa_func = testfunc;

left_bound = 0;
right_bound = 11;
raw_err_tol = 1e-1;

% left_root = vpasolve(vpa_func == 0, vpa_var, left_bound);
% right_root = vpasolve(vpa_func == 0, vpa_var, right_bound);
left_root = fzero(vpa_func, left_bound);
right_root = fzero(vpa_func, right_bound);
known_roots_vec = [left_root, right_root];
roots_vec = find_roots_binsearch(vpa_func, vpa_var, left_bound, right_bound, left_root, right_root, raw_err_tol, known_roots_vec);