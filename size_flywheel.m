function [D_outer, I_required, mass, energy_fluctuation] = size_flywheel(T_total, theta, params)
    % SIZE_FLYWHEEL Calculate required flywheel dimensions based on energy fluctuation
    %
    % Inputs:
    %   T_total - total torque array (N*m)
    %   theta - crank angle array (radians)
    %   params - engine parameters structure
    %
    % Outputs:
    %   D_outer - required outer diameter (m)
    %   I_required - required moment of inertia (kg*m^2)
    %   mass - flywheel mass (kg)
    %   energy_fluctuation - maximum energy fluctuation (J)
    
    % Extract parameters
    omega_avg = params.omega_avg;              % Average angular velocity (rad/s)
    Cs = params.flywheel.coeff_fluctuation;    % Coefficient of speed fluctuation
    w = params.flywheel.width;                 % Flywheel width (m)
    t = params.flywheel.thickness;             % Rim thickness (m)
    rho = params.flywheel.material_density;    % Material density (kg/m^3)
    
    % Calculate mean torque
    T_mean = mean(T_total);
    
    % Calculate energy fluctuation throughout the cycle
    % Energy = integral of (T - T_mean) * d(theta)
    
    % Method 1: Direct integration of torque deviation
    T_deviation = T_total - T_mean;  % Torque deviation from mean
    
    % Calculate cumulative energy variation
    energy_variation = zeros(size(theta));
    for i = 2:length(theta)
        dtheta = theta(i) - theta(i-1);
        % Energy change = torque * angular displacement
        energy_variation(i) = energy_variation(i-1) + ...
                              0.5 * (T_deviation(i) + T_deviation(i-1)) * dtheta;
    end
    
    % Maximum energy fluctuation
    E_max = max(energy_variation);
    E_min = min(energy_variation);
    energy_fluctuation = E_max - E_min;
    
    % Required moment of inertia
    % From the coefficient of fluctuation formula:
    % Cs = (omega_max - omega_min) / omega_avg
    % Energy fluctuation = I * omega_avg^2 * Cs
    
    I_required = energy_fluctuation / (Cs * omega_avg^2);
    
    % For an annular (hollow) cylinder flywheel:
    % I = (1/2) * m * (r_outer^2 + r_inner^2)
    % where r_inner = r_outer - t
    
    % Volume of annular cylinder: V = pi * w * (r_outer^2 - r_inner^2)
    % Mass: m = rho * V = rho * pi * w * (r_outer^2 - r_inner^2)
    
    % Substituting r_inner = r_outer - t:
    % V = pi * w * (r_outer^2 - (r_outer - t)^2)
    % V = pi * w * (2*r_outer*t - t^2)
    % m = rho * pi * w * (2*r_outer*t - t^2)
    
    % For the moment of inertia:
    % I = (1/2) * m * (r_outer^2 + (r_outer - t)^2)
    % I = (1/2) * rho * pi * w * (2*r_outer*t - t^2) * (r_outer^2 + (r_outer - t)^2)
    
    % Simplifying for thin rim (t << r_outer):
    % I ≈ m * r_outer^2 ≈ rho * pi * w * 2 * r_outer * t * r_outer^2
    % I ≈ 2 * pi * rho * w * t * r_outer^3
    
    % Solving for r_outer:
    % r_outer = (I_required / (2 * pi * rho * w * t))^(1/3)
    
    % More accurate calculation without thin-rim assumption
    % This requires solving a cubic equation, but we can iterate
    
    % Initial guess using thin-rim approximation
    r_outer_guess = (I_required / (2 * pi * rho * w * t))^(1/3);
    
    % Iterative refinement
    for iter = 1:10
        r_inner = r_outer_guess - t;
        if r_inner <= 0
            error('Flywheel thickness too large for required inertia');
        end
        
        % Calculate mass with current guess
        V = pi * w * (r_outer_guess^2 - r_inner^2);
        m = rho * V;
        
        % Calculate moment of inertia
        I_calc = 0.5 * m * (r_outer_guess^2 + r_inner^2);
        
        % Adjust radius based on error
        error_ratio = I_required / I_calc;
        r_outer_guess = r_outer_guess * error_ratio^(1/3);
        
        % Check convergence
        if abs(I_calc - I_required) / I_required < 0.001
            break;
        end
    end
    
    % Final values
    r_outer = r_outer_guess;
    D_outer = 2 * r_outer;  % Diameter
    r_inner = r_outer - t;
    
    % Calculate final mass
    V_flywheel = pi * w * (r_outer^2 - r_inner^2);
    mass = rho * V_flywheel;
    
    % Verify moment of inertia
    I_final = 0.5 * mass * (r_outer^2 + r_inner^2);
    
    % Calculate kinetic energy stored
    KE_stored = 0.5 * I_final * omega_avg^2;
    
    % Display flywheel design results
    fprintf('Flywheel Design:\n');
    fprintf('  Energy Fluctuation: %.2f J\n', energy_fluctuation);
    fprintf('  Required Moment of Inertia: %.4f kg*m^2\n', I_required);
    fprintf('  Outer Diameter: %.3f m (%.1f mm)\n', D_outer, D_outer*1000);
    fprintf('  Inner Diameter: %.3f m (%.1f mm)\n', 2*r_inner, 2*r_inner*1000);
    fprintf('  Mass: %.2f kg\n', mass);
    fprintf('  Kinetic Energy (at avg speed): %.2f kJ\n', KE_stored/1000);
    fprintf('  Material: Steel (density = %.0f kg/m^3)\n', rho);
    fprintf('  Width: %.3f m\n', w);
    fprintf('  Rim Thickness: %.3f m\n', t);
    fprintf('\n');
    
    % Validate design
    if D_outer > params.limits.max_flywheel_diameter
        warning('Flywheel diameter (%.2f m) exceeds maximum limit (%.2f m)', ...
                D_outer, params.limits.max_flywheel_diameter);
    end
    
    if D_outer < 0.1
        warning('Flywheel diameter seems too small (%.3f m) - check parameters', D_outer);
    end
    
    if mass > 1000
        warning('Flywheel mass exceeds 1000 kg - may be impractical');
    end
    
    % Check if moment of inertia requirement is met
    inertia_error = abs(I_final - I_required) / I_required * 100;
    if inertia_error > 1
        warning('Moment of inertia error: %.2f%%', inertia_error);
    end
end