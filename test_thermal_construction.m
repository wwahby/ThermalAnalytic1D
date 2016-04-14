% test constructor

close all
clear all

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

% alpha = [alpha_si, alpha_ox, alpha_si];
% k_actual = [k_si, k_ox, k_si];
% thickness_actual = [50, 5, 50] * 1e-6;
% pdens_cm2 = [0, 100, 0];

dTR = 1;

%%
cooling_configs = {[h_air, h_package], [h_water, h_package]};
tiers = [1, 2, 3];
%tier_scaling_funcs = {@(g,n) 1, @(g,n) g*(1+(1/g-1)./log(exp(1)+n-1))};
tier_scaling_funcs = {@(g,n) g*(1+(1/g-1)./log(exp(1)+n-1))};
pdens_cm2_base = 100;


num_cooling_configs = length(cooling_configs);
num_tiers = length(tiers);
num_scaling_funcs = length(tier_scaling_funcs);

theta_cell = cell(num_cooling_configs, num_tiers, num_scaling_funcs);
eta_max_mat = zeros(num_tiers);

alpha_layers = [alpha_si, alpha_ox];
k_layers = [k_si, k_ox];
thickness_layers = [50, 5]*1e-6;
pdens_layers = [0, pdens_cm2_base];

gg = 0.5;


for cind = 1:num_cooling_configs
    h_actual = cooling_configs{cind};
    
    for nind = 1:num_tiers
        n = tiers(nind);
        
        % Construct inputs
        alpha = zeros(1, 2*n);
        k_actual = zeros(1, 2*n);
        pdens_cm2_base_vec = zeros(1, 2*n);
        thickness_actual = zeros(1, 2*n);
        for ii = 1:n
            alpha(2*ii-1:2*ii) = alpha_layers(1:2);
            k_actual(2*ii-1:2*ii) = k_layers(1:2);
            thickness_actual(2*ii-1:2*ii) = thickness_layers(1:2);
            pdens_cm2_base_vec(2*ii-1:2*ii) = pdens_layers(1:2);
        end
        
        for scind = 1:num_scaling_funcs
            scaling_func = tier_scaling_funcs{scind};
            pdens_scaling_factor = scaling_func(gg, n);
            pdens_cm2 = pdens_cm2_base_vec * pdens_scaling_factor;
            [theta, eta_vec, ~] = construct_thermal_function_cartesian(alpha, k_actual, thickness_actual, pdens_cm2, h_actual);
            theta_cell{cind, nind, scind} = theta;
            eta_max_mat(nind) = max(eta_vec);
        end
    end
end


%%
Nx = 1e3;
Nt = 1e3;
tvec = logspace(-6,0,Nt);
dt_max_mat = zeros(num_cooling_configs, num_tiers, num_scaling_funcs, Nt);

for cind = 1:num_cooling_configs
    for nind = 1:num_tiers
        for scind = 1:num_scaling_funcs
            for tind = 1:Nt
                xvec = linspace(0, 0.99*eta_max_mat(nind), Nx);
                dt_max_mat(cind, nind, scind, tind) = max( theta_cell{cind, nind, scind}(xvec, tvec(tind)) );
            end
        end
    end
end

%%
colors = {'b', 'r'};
figure(1)
clf
hold on
for cind=1:num_cooling_configs
    for nind = 1:num_tiers
        for scind = 1:num_scaling_funcs
            dt_max_vec = squeeze(dt_max_mat(cind,nind,scind,:));
            plot(tvec, dt_max_vec, 'color', colors{cind});
        end
    end
end
xlabel('Time (s)')
ylabel('Max \DeltaT (K)')
set(gca,'xscale','log')
fixfigs(1,3,14,12)

figure(2)
clf
hold on
for cind=1:num_cooling_configs
    for nind = 1:num_tiers
        for scind = 1:num_scaling_funcs
            dt_max_vec = squeeze(dt_max_mat(cind,nind,scind,:));
            plot(tvec, dt_max_vec, 'color', colors{cind});
        end
    end
end
xlabel('Time (s)')
ylabel('Max \DeltaT (K)')
set(gca,'xscale','log')
set(gca,'yscale','log')
fixfigs(2,3,14,12)

figure(3)
clf
hold on
for cind=1:num_cooling_configs
    for nind = 1:num_tiers
        for scind = 1:num_scaling_funcs
            dt_max_vec = squeeze(dt_max_mat(cind,nind,scind,:));
            plot(tvec, dt_max_vec, 'color', colors{cind});
        end
    end
end
xlabel('Time (s)')
ylabel('Max \DeltaT (K)')
xlim([1e-3,2e-2])
ylim([0, 60])
set(gca,'xscale','log')
fixfigs(3,3,14,12)
