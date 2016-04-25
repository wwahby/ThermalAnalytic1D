% test constructor

close all
clear all

scalar_func = @(x) ones(1,length(x));
scalar_zero_func = @(x) zeros(1,length(x));


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

h_actual = [h_air, h_air];


% alpha = [alpha_si, alpha_ox, alpha_si, alpha_ox, alpha_si, alpha_ox, alpha_si, alpha_ox, alpha_si, alpha_ox];
% k_actual = [k_si, k_ox, k_si, k_ox, k_si, k_ox, k_si, k_ox, k_si, k_ox];
% thickness_actual = [50, 5, 50, 5, 50, 5, 50, 5, 50, 5] * 1e-6;
% pdens_cm2 = [0, 100, 0, 100, 0, 100, 0, 100, 0, 100];

% too big -- det conversion fails because of excessive { [ ( nesting in
% str2func. "Error: Nesting of {, [, and ( cannot exceed a depth of 32."
% alpha = [alpha_si, alpha_ox, alpha_si, alpha_ox, alpha_si, alpha_ox, alpha_si, alpha_ox, alpha_si];
% k_actual = [k_si, k_ox, k_si, k_ox, k_si, k_ox, k_si, k_ox, k_si];
% thickness_actual = [50, 5, 50, 5, 50, 5, 50, 5, 50] * 1e-6;
% pdens_cm2 = [0, 100, 0, 100, 0, 100, 0, 100, 0];

% alpha = [alpha_si, alpha_ox, alpha_si, alpha_ox, alpha_si, alpha_ox, alpha_si, alpha_ox];
% k_actual = [k_si, k_ox, k_si, k_ox, k_si, k_ox, k_si, k_ox];
% thickness_actual = [50, 5, 50, 5, 50, 5, 50, 5] * 1e-6;
% pdens_cm2 = [0, 100, 0, 100, 0, 100, 0, 100];

% alpha = [alpha_si, alpha_ox, alpha_si, alpha_ox, alpha_si, alpha_ox, alpha_si];
% k_actual = [k_si, k_ox, k_si, k_ox, k_si, k_ox, k_si];
% thickness_actual = [50, 5, 50, 5, 50, 5, 50] * 1e-6;
% pdens_cm2 = [0, 100, 0, 100, 0, 100, 0];

% alpha = [alpha_si, alpha_ox, alpha_si, alpha_ox, alpha_si, alpha_ox];
% k_actual = [k_si, k_ox, k_si, k_ox, k_si, k_ox];
% thickness_actual = [50, 5, 50, 5, 50, 5] * 1e-6;
% pdens_cm2 = [0, 100, 0, 100, 0, 100];

% alpha = [alpha_si, alpha_ox, alpha_si, alpha_ox, alpha_si];
% k_actual = [k_si, k_ox, k_si, k_ox, k_si];
% thickness_actual = [50, 5, 50, 5, 50] * 1e-6;
% pdens_cm2 = [0, 100, 0, 100, 0];

% alpha = [alpha_si, alpha_ox, alpha_si, alpha_ox];
% k_actual = [k_si, k_ox, k_si, k_ox];
% thickness_actual = [50, 5, 50, 5] * 1e-6;
% pdens_cm2 = [0, 100, 0, 100];
% power_functions = { @(x) 0, @(x) x, @(x) 0, @(x) x};

alpha = [alpha_si, alpha_ox, alpha_si];
k_actual = [k_si, k_ox, k_si];
thickness_actual = [50, 5, 50] * 1e-6;
pdens_cm2 = [0, 100, 0];
power_functions = {scalar_zero_func, scalar_func, scalar_zero_func};


% alpha = [alpha_si];
% k_actual = [k_si];
% thickness_actual = [20] * 1e-6;
% pdens_cm2 = [20];
% power_functions = { @(x) 1 + 100*(x>15e-6)};



dTR = 1;

dimensions = 1;
%[theta, eta_vec, ~, normalized_power_functions] = construct_thermal_function_cartesian(alpha, k_actual, thickness_actual, pdens_cm2, power_functions, h_actual);
[theta, eta_vec, ~, normalized_power_functions] = construct_thermal_function(alpha, k_actual, thickness_actual, pdens_cm2, power_functions, h_actual,dimensions);

%%
xvec = linspace(0, 0.999*max(eta_vec), 1e3);

figure(1)
clf
hold on
plot(xvec, theta(xvec, 1e-6))
plot(xvec, theta(xvec, 1e-5))
plot(xvec, theta(xvec, 1e-4))
xlim([0 max(eta_vec)])
xlabel('Normalized Position (-)')
ylabel('\DeltaT (K)')
fixfigs(1,3,14,12)

figure(2)
clf
hold on
plot(xvec, theta(xvec, 1e-3))
plot(xvec, theta(xvec, 1e-2))
plot(xvec, theta(xvec, 1e-1))
xlim([0 max(eta_vec)])
xlabel('Normalized Position (-)')
ylabel('\DeltaT (K)')
fixfigs(2,3,14,12)