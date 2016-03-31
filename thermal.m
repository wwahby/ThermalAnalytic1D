%% User input

k_si = 130;
k_ox = 1.38;
k_cu = 385;

alpha_si = k_si/2329/700;
alpha_ox = k_ox/2203/703;
alpha_cu = 1.11e-4;

h_air = 1.8e4;
h_water = 4.6e4;
h_package = 5;

alpha = [alpha_si, alpha_ox, alpha_si, alpha_ox, alpha_si];
k_actual = [k_si, k_ox, k_si, k_ox, k_si];
thickness_actual = [50, 5, 50, 5, 50] * 1e-6;
pdens_cm2 = [0, 100, 0, 100, 0];

% alpha = [alpha_si, alpha_ox, alpha_si];
% k_actual = [k_si, k_ox, k_si];
% thickness_actual = [50, 5, 50] * 1e-6;
% pdens_cm2 = [0, 100, 0];

h_actual = [h_air, h_water];

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

%% Symbolic stuff
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

%% Construct Matrix

A = N1;
for mind = 1:num_layers-1
        A = [A; P(mind, eta_vec, a, b_vec) ; Q(mind, eta_vec, a, b_vec, k_vec)];
end
A = [A; PM];
%%
detA = det(A);

fprintf('Finding roots...');
roots_start = cputime;

num_guesses = 3e1;
min_guess = 0;
max_guess = 5;
roots_vec = -1*ones(1,num_guesses);
init_guess_vec = linspace(min_guess, max_guess, num_guesses);
for gind = 1:length(init_guess_vec)
    if mod(gind, 5) == 0
        fprintf('%d\t',gind)
    end
    guess = init_guess_vec(gind);
    roots_vec(gind) = vpasolve(vpa(detA) == 0, a, guess);
    %roots_vec(gind) = solve(detA == 0, a, guess);
end
fprintf('\n')

roots_vec = roots_vec(roots_vec > 0);
roots_vec = unique(roots_vec);

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
wfunc = (x^p*wvec)';

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
xmin_lim = [0 eta_vec(1:end-1)];
xmax_lim = eta_vec;

denom_mat = sym(zeros(num_layers,length(roots_vec)) );
gnom_mat = sym(zeros(num_layers,length(roots_vec)) );
enom_mat = sym(zeros(num_layers,length(roots_vec)) );

int_time_start = cputime;
% Numeric integration -- faster than symbolic
for iind=1:num_layers
    for rind = 1:length(roots_vec)
        denom_mat(iind,rind) = integral( matlabFunction(WR2mat(iind,rind)), xmin_lim(iind), xmax_lim(iind));
        enom_mat(iind,rind) = integral( matlabFunction(WRmat(iind,rind)), xmin_lim(iind), xmax_lim(iind));
        if g_vec(iind) > 0
            gnom_mat(iind,rind) = integral( matlabFunction(gWRmat(iind,rind)), xmin_lim(iind), xmax_lim(iind));
        else
            gnom_mat(iind,rind) = 0;
        end
            
    end
end

% Symbolic integration -- ~7X slower than numerical integration
% for iind=1:num_layers
%     denom_mat(iind,:) = int( WR2mat(iind,:), x, [xmin_lim(iind), xmax_lim(iind)]);
%     gnom_mat(iind,:) = int( gWRmat(iind,:), x, [xmin_lim(iind), xmax_lim(iind)]);
%     enom_mat(iind,:) = int( WRmat(iind,:), x, [xmin_lim(iind), xmax_lim(iind)]);
% end
int_time_stop = cputime;
fprintf('\t(%.3g s)\n', int_time_stop - int_time_start);
%%

denom_n = sum(denom_mat,1);
gnum_n = sum(gnom_mat, 1);
enum_n = sum(enom_mat, 1);

gnom_n = gnum_n ./ denom_n;
enom_n = enum_n ./ denom_n;

%%
fprintf('Constructing final functions...\n');
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

%%
theta = @(x,t) theta5(alpha(1)*t/x_actual(1)^2, x);

xvec = linspace(0, 0.99*max(eta_vec), 1e3);


%% 

theta(xvec, 1e3);

figure(1)
clf
hold on
plot(xvec, theta(xvec, 1e-6))
plot(xvec, theta(xvec, 1e-5))
plot(xvec, theta(xvec, 1e-4))
xlim([0 max(eta_vec)])
fixfigs(1,3,14,12)

figure(2)
clf
hold on
plot(xvec, theta(xvec, 1e-3))
plot(xvec, theta(xvec, 1e-2))
plot(xvec, theta(xvec, 1e-1))
xlim([0 max(eta_vec)])
fixfigs(2,3,14,12)