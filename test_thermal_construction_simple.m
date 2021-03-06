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

% Materials: % Si, Ox, Cu
material_alphas = [alpha_si, alpha_ox, alpha_cu];
material_thermal_conds = [k_si, k_ox, k_cu];

h_air = 1.8e4;
h_water = 4.6e4;
h_package = 5;

h_actual = [h_package, h_air];

% materials = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2];
% thickness_actual = [50, 5, 50, 5, 50, 5, 50, 5, 50, 5] * 1e-6;
% pdens_cm2 = [0, 100, 0, 100, 0, 100, 0, 100, 0, 100];
% dimensions = 1;

% too big for normal conversion -- det conversion fails because of excessive { [ ( nesting in
% str2func. "Error: Nesting of {, [, and ( cannot exceed a depth of 32."
% materials = [1, 2, 1, 2, 1, 2, 1, 2, 1];
% thickness_actual = [50, 5, 50, 5, 50, 5, 50, 5, 50] * 1e-6;
% pdens_cm2 = [0, 100, 0, 100, 0, 100, 0, 100, 0];
% dimensions = 1;

% materials = [1, 2, 1, 2, 1, 2, 1, 2];
% thickness_actual = [50, 5, 50, 5, 50, 5, 50, 5] * 1e-6;
% pdens_cm2 = [0, 100, 0, 100, 0, 100, 0, 100];
% dimensions = 1;

% materials = [1, 2, 1, 2, 1, 2, 1];
% thickness_actual = [50, 5, 50, 5, 50, 5, 50] * 1e-6;
% pdens_cm2 = [0, 100, 0, 100, 0, 100, 0];
% dimensions = 1;

% materials = [1, 2, 1, 2, 1, 2];
% thickness_actual = [50, 5, 50, 5, 50, 5] * 1e-6;
% pdens_cm2 = [0, 100, 0, 100, 0, 100];
% dimensions = 1;

% materials = [1, 2, 1, 2, 1];
% thickness_actual = [50, 5, 50, 5, 50] * 1e-6;
% pdens_cm2 = [0, 100, 0, 100, 0];
% dimensions = 1;

% materials = [1, 2, 1, 2];
% thickness_actual = [50, 5, 50, 5] * 1e-6;
% pdens_cm2 = [0, 100, 0, 100];
% dimensions = 1;
 
% materials = [1, 2, 1];
% thickness_actual = [50, 5, 50] * 1e-6;
% pdens_cm2 = [0, 100, 0];
% dimensions = 1;

% materials = [3, 2, 1]; % cu, ox, si
% thickness_actual = [8, 0.1, 50] * 1e-6;
% pdens_cm2 = [100, 0, 0];
% dimensions = 2;


materials = [1, 2, 3];
thickness_actual = [50, 5, 50] * 1e-6;
pdens_cm2 = [0, 100, 0];
dimensions = 1;

materials = [2, 1, 2, 1];
thickness_actual = [5, 1, 5, 50] * 1e-6;
pdens_cm2 = [60, 0, 60, 0];
dimensions = 1;

% materials = [2, 1, 2, 1, 2, 1];
% thickness_actual = [5, 1, 5, 1, 5, 50] * 1e-6;
% pdens_cm2 = [75, 0, 75, 0, 75, 0];
% dimensions = 1;

%% Setup
power_functions = cell(1,length(pdens_cm2) );
for pind = 1:length(pdens_cm2)
    if (pdens_cm2(pind) > 0)
        power_functions{pind} = scalar_func;
    else
        power_functions{pind} = scalar_zero_func;
    end
end

alpha = material_alphas(materials);
k_actual = material_thermal_conds(materials);


[theta, eta_vec, ~, normalized_power_functions] = construct_thermal_function(alpha, k_actual, thickness_actual, pdens_cm2, power_functions, h_actual, dimensions);

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
plot(xvec, theta(xvec, 2e-2))
plot(xvec, theta(xvec, 4e-2))
plot(xvec, theta(xvec, 1e-1))
xlim([0 max(eta_vec)])
xlabel('Normalized Position (-)')
ylabel('\DeltaT (K)')
fixfigs(2,3,14,12)

tvec = logspace(-6,-1,5e1);
tmax = zeros(1, length(tvec));
for tind = 1:length(tvec)
    tmax(tind) = max(theta(xvec, tvec(tind)) );
end

figure(3)
clf
hold on
plot(tvec, tmax)
set(gca, 'yscale','log')
set(gca, 'xscale','log')
xlabel('time (s)')
ylabel('\DeltaT (K)')
grid on
fixfigs(3,3,14,12)
    