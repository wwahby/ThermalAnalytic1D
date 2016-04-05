function roots_vec = find_all_roots_in_bounds_fzero(func, bound_vec, Npts)
%func = @(x) sin(x);
% bound_vec = [0, 10];
% Npts = 3e2;

dx = (bound_vec(2) - bound_vec(1))/Npts;

pt_vec = bound_vec(1):dx:bound_vec(2);
val_vec = func(pt_vec);
sign_vec = sign(val_vec);

% actual_zeros = find(sign_vec == 0);
sign_change_inds = find( diff(sign_vec) );
bound_vecs = [sign_change_inds ; (sign_change_inds+1)]';
roots_vec = [];

for iind = 1:length(bound_vecs)
    if mod(iind,5) == 0;
        fprintf(' %d', iind);
    end
    bounds = pt_vec(bound_vecs(iind, :));
    new_root = fzero(func, bounds);
    roots_vec = [roots_vec, new_root];
end
