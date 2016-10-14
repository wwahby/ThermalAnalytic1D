% test constructor

close all
clear all

%%
scalar_func = @(x) ones(1,length(x));
scalar_zero_func = @(x) zeros(1,length(x));

tvec = logspace(-3,2,5e2);
tlimit = 70;

%% Constants and material parameters
k_si = 130;
k_ox = 1.38;
k_cu = 385;

alpha_si = k_si/2329/700;
alpha_ox = k_ox/2203/703;
alpha_cu = 1.11e-4;

% Materials: % Si, Ox, Cu
material_alphas = [alpha_si, alpha_ox, alpha_cu];
material_thermal_conds = [k_si, k_ox, k_cu];

h_passive = 1e3;
h_air = 1.8e4;
h_water = 4.6e4;
h_package = 5;

h_actual = [h_package, h_air];

%% 2D
hs_thickness_um = 4e2;

materials = [2, 1, 3];
thickness_actual = [5, 200, hs_thickness_um] * 1e-6;
pdens_cm2 = [100, 0, 0];
dimensions = 1;

alpha = material_alphas(materials);
k_actual = material_thermal_conds(materials);
tmax_2d = zeros(1, length(tvec));

[theta, eta_vec_2d, A, normalized_power_functions] = construct_thermal_function_uniform(alpha, k_actual, thickness_actual, pdens_cm2, h_actual, dimensions);
theta_2d = theta;


%% 3DIC

materials = [2, 1, 2, 2, 1, 3];
thickness_actual = [5, 200, 10, 5, 200, hs_thickness_um] * 1e-6;
pdens_cm2 = [100, 0, 0, 10, 0, 0];
dimensions = 1;

alpha = material_alphas(materials);
k_actual = material_thermal_conds(materials);
tmax_3d = zeros(1, length(tvec));

[theta, eta_vec_3d, A, normalized_power_functions] = construct_thermal_function_uniform(alpha, k_actual, thickness_actual, pdens_cm2, h_actual, dimensions);
theta_3d = theta;


%% 3DIC B2B

materials = [2, 1, 2, 1, 2, 3];
thickness_actual = [5, 200, 10, 200, 5, hs_thickness_um] * 1e-6;
pdens_cm2 = [100, 0, 0, 0, 10, 0];
dimensions = 1;

alpha = material_alphas(materials);
k_actual = material_thermal_conds(materials);
tmax_3db = zeros(1, length(tvec));

[theta, eta_vec_3db, A, normalized_power_functions] = construct_thermal_function_uniform(alpha, k_actual, thickness_actual, pdens_cm2, h_actual, dimensions);
theta_3db = theta;

%% Times to evaluate and max temp


% Get tmax and time to exceed limit

for time_ind = 1:length(tvec)
    xvec = linspace(0, 0.999*max(eta_vec_2d), 1e3);
    tmax_2d(time_ind) = max(theta_2d(xvec, tvec(time_ind)) );
end
time_ind = find(tmax_2d < tlimit, 1,'last');
time_to_limit_2d = tvec(time_ind);

for time_ind = 1:length(tvec)
    xvec = linspace(0, 0.999*max(eta_vec_3d), 1e3);
    tmax_3d(time_ind) = max(theta_3d(xvec, tvec(time_ind)) );
end
time_ind = find(tmax_3d < tlimit, 1,'last');
time_to_limit_3d = tvec(time_ind);

for time_ind = 1:length(tvec)
    xvec = linspace(0, 0.999*max(eta_vec_3db), 1e3);
    tmax_3db(time_ind) = max(theta_3db(xvec, tvec(time_ind)) );
end
time_ind = find(tmax_3db < tlimit, 1,'last');
time_to_limit_3db = tvec(time_ind);

%% Plots
figure(1)
clf
hold on
plot(tvec, tmax_2d, 'k')
plot(tvec, tmax_3d, 'b')
plot(tvec, tmax_3db, 'r')
grid on
xlabel('time (s)')
ylabel('\DeltaT (K)')
set(gca,'xscale','lin')
xlim([0 6])
fixfigs(1,2,14,12)

figure(2)
clf
hold on
plot(tvec, tmax_2d, 'k')
plot(tvec, tmax_3d, 'b')
plot(tvec, tmax_3db, 'r')
grid on
xlabel('time (s)')
ylabel('\DeltaT (K)')
set(gca,'xscale','log')
xlim([1e-3, 1e1])
fixfigs(2,2,14,12)