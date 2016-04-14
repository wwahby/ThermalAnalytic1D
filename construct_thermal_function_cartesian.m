function [theta, eta_vec, A] = construct_thermal_function_cartesian(alpha, k, thicknesses, pdens_cm2, h)
time_start = cputime;

thickness_actual = thicknesses;
k_actual = k;
h_actual = h;
dTR = 1;

%% Automatic Setup

x_actual = cumsum(thickness_actual);

dx = [ x_actual(1) diff(x_actual)];

% Construct useful constants
pdens_m2 = pdens_cm2 * 1e4;
pdens_m3 = pdens_m2 ./ dx;
hM = h_actual(end) * x_actual(1)/k_actual(end);
h1 = h_actual(1) * x_actual(1)/k_actual(1);

num_layers = length(alpha);

% Construct important vectors
k_vec = k_actual(2:end)./k_actual(1:end-1);
eta_vec = x_actual./x_actual(1);
b_vec = alpha ./ alpha(1);
g_vec = b_vec .* x_actual(1)^2 .* pdens_m3./ k_actual  / dTR;

%% Symbolic stuff - -Construct basis functions

fprintf('Constructing Matrix...');
mat_start = cputime;

syms a x b;

phis = cos(a*x/sqrt(b));
psis = sin(a*x/sqrt(b));

dphis = diff(phis,x);
dpsis = diff(psis,x);

phi = @(xx, aa, bb) subs(phis, [x, a, b], [xx, aa, bb]);
psi = @(xx, aa, bb) subs(psis, [x, a, b], [xx, aa, bb]);
dphi = @(xx, aa, bb) subs(dphis, [x, a, b], [xx, aa, bb]);
dpsi = @(xx, aa, bb) subs(dpsis, [x, a, b], [xx, aa, bb]);

N1 = [ -dphi(0, a, b_vec(1))/h1 + phi(0, a, b_vec(1)), -dpsi(0, a, b_vec(1))/h1 + psi(0, a, b_vec(1)), zeros(1,2*(num_layers-1)) ];
P = @(n, eta, a, b_vec) [ zeros(1, 2*(n-1)), phi(eta(n), a, b_vec(n)), psi(eta(n), a, b_vec(n)), -phi( eta(n), a, b_vec(n+1) ), -psi( eta(n), a, b_vec(n+1) ), zeros(1, 2*(num_layers-n-1)) ];
Q = @(n, eta, a, b_vec, k_vec) [ zeros(1, 2*(n-1)), dphi(eta(n), a, b_vec(n)), dpsi(eta(n), a, b_vec(n)), -k_vec(n)*dphi( eta(n), a, b_vec(n+1) ), -k_vec(n)*dpsi( eta(n), a, b_vec(n+1) ), zeros(1, 2*(num_layers-n-1)) ];
PM = [ zeros(1, 2*(num_layers-1)), dphi(eta_vec(end), a, b_vec(end))/hM + phi(eta_vec(end), a, b_vec(end)), dpsi( eta_vec(end), a, b_vec(end))/hM + psi( eta_vec(end), a, b_vec(end))];

%% Build Matrix
A = N1;
for mind = 1:num_layers-1
        A = [A; P(mind, eta_vec, a, b_vec) ; Q(mind, eta_vec, a, b_vec, k_vec)];
end
A = [A; PM];

mat_stop = cputime;
fprintf('\t(%.3g s)\n', mat_stop - mat_start);


%% Take determinant -- alt det -- 50X faster
fprintf('Taking determinant...');
det_start = cputime;
[L, U, P] = lu(vpa(A));
detA = det(P^-1)*prod(diag(L))*prod(diag(U)); % using the fact that det(AB) = det(A)*det(B) and the fact that det(triangular matrix) = prod(diagonal entries)
det_stop = cputime;
fprintf('\t(%.3g s)\n', det_stop - det_start);

%% Third alternate determinant conversion method
fprintf('Converting determinant (sym->str)...');
det3_start = cputime;
detstr = ['@(a) ', char(detA)]; % incredibly hacky way to quickly convert symbolic determinant into matlab function handle. matlabFunction takes FOREVER for long symbolics.
detstr_a = strrep( detstr, '*', '.*');
detstr_b = strrep(detstr_a, '^', '.^');
detstr_c = strrep(detstr_b, '/', './');

%detfunc = eval(detstr_c);
trigger_depth = 20;
[detfunc, ~, ~] = condense_functions(detstr_c, trigger_depth);
det3_stop = cputime;
fprintf('\t(%.3g s)\n', det3_stop - det3_start);

%% Alternate root method
fprintf('Finding roots...');
roots_start = cputime;
bound_vec = [1e-30,10]; % have to use something greater than zero because the LU determinant method gives us something with a NaN at precisely zero, but matches everywhere else (and is 50X faster)
Npts = 3e2;
roots_vec = find_all_roots_in_bounds_fzero(detfunc, bound_vec, Npts);
roots_vec = unique( roots_vec(roots_vec > 0) );

roots_stop = cputime;
fprintf('\t(%.3g s)\n', roots_stop - roots_start);

%%
fprintf('Finding coefficients...');
coeff_start = cputime;
null_list = zeros(2*num_layers, length(roots_vec));
norm_list = zeros(2*num_layers, length(roots_vec));
for rind = 1:length(roots_vec)
    root = roots_vec(rind);
    null_list(:,rind) = null( double(vpa(subs(A,a,root))) );
    
    % Normalize nullspace basis
    norm_list(:,rind) = null_list(:,rind)/null_list(1,rind);
end

C_mat = norm_list(1:2:end,:);
D_mat = norm_list(2:2:end,:);

coeff_stop = cputime;
fprintf('\t(%.3g s)\n', coeff_stop - coeff_start);

%% Construct Phi/Psi vectors
fprintf('Constructing basis vectors...');
basis_start = cputime;
phi_vec_base = subs(phi(x,a,b),b,b_vec)';
psi_vec_base = subs(psi(x,a,b),b,b_vec)';

phi_func_mat = sym(zeros(num_layers,length(roots_vec)) );
psi_func_mat = sym(zeros(num_layers,length(roots_vec)) );

for iind=1:num_layers
   phi_func_mat(iind,:) = C_mat(iind,:) .*phi_vec_base(iind);
   psi_func_mat(iind,:) = D_mat(iind,:) .*psi_vec_base(iind);
end
R_func_mat = phi_func_mat + psi_func_mat;

R_vec = R_func_mat;
for rind = 1:length(roots_vec)
    R_vec(:,rind) = subs(R_vec(:,rind), a, roots_vec(rind));
end

basis_stop = cputime;
fprintf('\t(%.3g s)\n', basis_stop - basis_start);

%%
fprintf('Constructing weighting vectors...\t');
weighting_start = cputime;
Fvec = cumprod([1, k_vec]);
wvec = Fvec./b_vec;

p = 0;
wfunc = (x.^p*wvec)';

WRmat = sym(zeros(num_layers,length(roots_vec)) );
WR2mat = sym(zeros(num_layers,length(roots_vec)) );
gWRmat = sym(zeros(num_layers,length(roots_vec)) );

for iind = 1:num_layers
    WRmat(iind,:) = wfunc(iind).*R_vec(iind,:);
    WR2mat(iind,:) = wfunc(iind).*(R_vec(iind,:).^2);
    gWRmat(iind,:) = wfunc(iind).*g_vec(iind).*R_vec(iind,:);
end

weighting_stop = cputime;
fprintf('\t(%.3g s)\n', weighting_stop - weighting_start);

%% 
fprintf('Integrating weighting vectors...');
int_time_start = cputime;
xmin_lim = [0 eta_vec(1:end-1)];
xmax_lim = eta_vec;

R_cell = cell(num_layers, length(roots_vec));
w_cell = cell(1,num_layers);

denom_mat = zeros(num_layers,length(roots_vec));
gnom_mat = zeros(num_layers,length(roots_vec)) ;
enom_mat = zeros(num_layers,length(roots_vec)) ;

for iind = 1:num_layers
    w_cell{iind} = @(x) x.^p * wvec(iind);
    for rind = 1:length(roots_vec)
        R_cell{iind, rind} = str2func( ['@(x) ', char(R_vec(iind, rind)) ] );
        
        denom_mat(iind,rind) = integral( @(x) w_cell{iind}(x).*(R_cell{iind, rind}(x).^2), xmin_lim(iind), xmax_lim(iind) );
        enom_mat(iind,rind) = integral( @(x) w_cell{iind}(x).*(R_cell{iind, rind}(x)), xmin_lim(iind), xmax_lim(iind) );
        if g_vec(iind) > 0
            gnom_mat(iind,rind) = integral( @(x) g_vec(iind).*w_cell{iind}(x).*(R_cell{iind, rind}(x)), xmin_lim(iind), xmax_lim(iind) );
        else
            gnom_mat(iind,rind) = 0;
        end
    end
end


int_time_stop = cputime;
fprintf('\t(%.3g s)\n', int_time_stop - int_time_start);


denom_n = sum(denom_mat,1);
gnum_n = sum(gnom_mat, 1);
enum_n = sum(enom_mat, 1);

gnom_n = gnum_n ./ denom_n;
enom_n = enum_n ./ denom_n;

%%
fprintf('Constructing final functions...');
func_time_start = cputime;

syms Fo;
omega = 0;
AN = enom_n .* ( omega^2 ./ (omega^2 + roots_vec.^4) - 1) - gnom_n./roots_vec.^2;
ANsin = enom_n .* (omega * roots_vec.^2 ./ (omega^2 + roots_vec.^4)) - gnom_n./roots_vec.^2;
ANsinterm = ANsin.*exp(-roots_vec.^2 * Fo);
ENsinterm = enom_n .* ((omega^2*sin(omega*Fo) + roots_vec.^2*omega*cos(omega*Fo))/(omega^2 + roots_vec.^4) );
GNterm = gnom_n./roots_vec.^2;

mod_sin_term = ANsinterm - ENsinterm + GNterm;

thetafunc2sin = sym(zeros(num_layers,length(roots_vec)) );
limfunc = sym(zeros(num_layers,1));
sympref('HeavisideAtOrigin', 1); % need this because otherwise heaviside(0) = 0.5
for iind = 1:num_layers
    thetafunc2sin(iind, :) = R_vec(iind,:) .* mod_sin_term;
    limfunc(iind) = heaviside(x-xmin_lim(iind)) - heaviside(x-xmax_lim(iind));
end

thetafunc3sin = sin(omega*Fo) + sum(thetafunc2sin, 2);
theta4 = sum(thetafunc3sin .* limfunc);

theta5 = matlabFunction(theta4);

theta = @(x,t) theta5(alpha(1)*t/x_actual(1)^2, x);

func_time_stop = cputime;
fprintf('\t(%.3g s)\n', func_time_stop - func_time_start);

time_stop = cputime;
fprintf('Total elapsed time: %.4g s (%.4g m)\n\n', time_stop - time_start, (time_stop - time_start)/60)