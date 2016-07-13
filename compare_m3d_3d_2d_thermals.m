% test constructor

close all
clear all

scalar_func = @(x) ones(1,length(x));
scalar_zero_func = @(x) zeros(1,length(x));


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

tvec = logspace(-6,-1,5e2);
tlimit = 70;

h_actual = [h_package, h_air];

%% 2D

materials = [2, 1];
thickness_actual = [5, 50] * 1e-6;
pdens_cm2 = [100, 0];
dimensions = 1;

alpha = material_alphas(materials);
k_actual = material_thermal_conds(materials);
tmax_2d = zeros(1, length(tvec));

[theta, eta_vec, A, normalized_power_functions] = construct_thermal_function_uniform(alpha, k_actual, thickness_actual, pdens_cm2, h_actual, dimensions);
theta_2d = theta;

for time_ind = 1:length(tvec)
    xvec = linspace(0, 0.999*max(eta_vec), 1e3);
    tmax_2d(time_ind) = max(theta(xvec, tvec(time_ind)) );
end
time_ind = find(tmax_2d < tlimit, 1,'last');
time_to_limit_2d = tvec(time_ind);

%% 3DIC

tiers = [2, 3, 4];
%pdens_cm2_logic = [65, 75, 85];
pdens_cm2_logic = (100 + tiers*300/8)./tiers;
thickness_ox = 5e-6;
thickness_si = 50e-6;

alpha = material_alphas(materials);
k_actual = material_thermal_conds(materials);

theta_3dic = cell(1,length(tiers));
time_to_limit_3dic = zeros(1, length(tiers));
tmax_3dic = zeros(length(tiers), length(tvec));

for tind = 1:length(tiers)
    pdens_logic = pdens_cm2_logic(tind);
    
    pdens_cm2 = zeros(1,2*tiers(tind));
    thickness_actual = zeros(1,2*tiers(tind));
    materials = zeros(1, 2*tiers(tind));
    
    pdens_cm2(1:2:end) = pdens_cm2_logic(tind);
    thickness_actual(1:2:end) = thickness_ox;
    thickness_actual(2:2:end) = thickness_si;
    materials(1:2:end) = 2;
    materials(2:2:end) = 1;
    
    alpha = material_alphas(materials);
    k_actual = material_thermal_conds(materials);
    
    [theta, eta_vec, A, normalized_power_functions] = construct_thermal_function_uniform(alpha, k_actual, thickness_actual, pdens_cm2, h_actual, dimensions);
    theta_3dic{tind} = theta;
    
    for time_ind = 1:length(tvec)
        xvec = linspace(0, 0.999*max(eta_vec), 1e3);
        tmax_3dic(tind, time_ind) = max(theta(xvec, tvec(time_ind)) );
    end
    time_last_ind = find(tmax_3dic(tind,:) < tlimit, 1,'last');
    time_to_limit_3dic(tind) = tvec(time_last_ind);
end


%% M3DIC

tiers = [2, 3, 4];
pdens_cm2_logic = (100 + tiers*300/8)./tiers;
thickness_ox = 5e-6;
thickness_si = 1e-6;
thickness_si_bulk = 50e-6;



theta_m3dic = cell(1,length(tiers));
time_to_limit_m3dic = zeros(1, length(tiers));
tmax_m3dic = zeros(length(tiers), length(tvec));

for tind = 1:length(tiers)
    pdens_logic = pdens_cm2_logic(tind);
    
    pdens_cm2 = zeros(1, 2*tiers(tind));
    thickness_actual = zeros(1, 2*tiers(tind));
    materials = zeros(1, 2*tiers(tind));
    
    pdens_cm2(1:2:end) = pdens_cm2_logic(tind);
    thickness_actual(1:2:end) = thickness_ox;
    thickness_actual(2:2:end) = thickness_si;
    thickness_actual(end) = thickness_si_bulk;
    materials(1:2:end) = 2;
    materials(2:2:end) = 1;
    
    alpha = material_alphas(materials);
    k_actual = material_thermal_conds(materials);
    
    [theta, eta_vec, A, normalized_power_functions] = construct_thermal_function_uniform(alpha, k_actual, thickness_actual, pdens_cm2, h_actual, dimensions);
    theta_m3dic{tind} = theta;
    
    for time_ind = 1:length(tvec)
        xvec = linspace(0, 0.999*max(eta_vec), 1e3);
        tmax_m3dic(tind, time_ind) = max(theta(xvec, tvec(time_ind)) );
    end
    time_last_ind = find(tmax_m3dic(tind,:) < tlimit, 1,'last');
    time_to_limit_m3dic(tind) = tvec(time_last_ind);
end 

%%

colors = {'r', 'b', 'g'};
linestyles = {'-', ':'};


figure(3)
clf
hold on
plot(tvec, tmax_2d, 'k-')
for tind = 1:length(tiers)
    plot(tvec, tmax_3dic(tind, :), 'color', colors{tind}, 'linestyle', '-')
    plot(tvec, tmax_m3dic(tind, :), 'color', colors{tind}, 'linestyle', '--')
end
set(gca, 'yscale','lin')
set(gca, 'xscale','log')
ylim([0 180])
ylim([0 70])
xlim([1e-3 1e-1])
xlabel('time (s)')
ylabel('\DeltaT (K)')
grid on
fixfigs(3,3,14,12)


time_to_limit_3d = [time_to_limit_2d time_to_limit_3dic];
time_to_limit_m3d = [time_to_limit_2d time_to_limit_m3dic];

figure(5)
clf
hold on
plot(time_to_limit_3d, 'b')
plot(time_to_limit_m3d, 'r')
xlabel('Logic tiers')
ylabel('Time to exceed \DeltaT_{max} (s)')
set(gca,'xtick',1:4)
set(gca,'yscale','log')
grid on
fixfigs(5,3,14,12);