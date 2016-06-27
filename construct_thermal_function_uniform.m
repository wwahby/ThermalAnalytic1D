function [theta, eta_vec, A, normalized_power_functions] = construct_thermal_function_uniform(alpha, k, thicknesses, pdens_cm2, h, dimensions)

scalar_func = @(x) ones(1,length(x));
scalar_zero_func = @(x) zeros(1,length(x));

power_functions = cell(1,length(pdens_cm2) );

for pind = 1:length(pdens_cm2)
    if (pdens_cm2(pind) > 0)
        power_functions{pind} = scalar_func;
    else
        power_functions{pind} = scalar_zero_func;
    end
end

[theta, eta_vec, A, normalized_power_functions] = construct_thermal_function(alpha, k, thicknesses, pdens_cm2, power_functions, h, dimensions);

