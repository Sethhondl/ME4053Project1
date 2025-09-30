function [omega, alpha, rpm, Cs_actual] = simulate_dynamics(T_total, theta, I_flywheel, params)
    % SIMULATE_DYNAMICS Simulate angular velocity variation with flywheel
    %
    % Inputs:
    %   T_total - total torque array (N*m)
    %   theta - crank angle array (radians)
    %   I_flywheel - flywheel moment of inertia (kg*m^2)
    %   params - engine parameters structure
    %
    % Outputs:
    %   omega - instantaneous angular velocity (rad/s)
    %   alpha - angular acceleration (rad/s^2)
    %   rpm - instantaneous speed in RPM
    %   Cs_actual - actual coefficient of speed fluctuation achieved
    
    % Extract parameters
    omega_avg = params.averageAngularVelocity;  % Average angular velocity (rad/s)
    T_mean = mean(T_total);        % Mean torque
    
    % Assume constant load torque equal to mean engine torque (steady state)
    T_load = T_mean;
    
    % Net torque (engine torque - load torque)
    T_net = T_total - T_load;
    
    % Angular acceleration: alpha = T_net / I
    alpha = T_net / I_flywheel;
    
    % Integrate acceleration to get velocity variation
    % Start with average angular velocity
    omega = zeros(size(theta));
    omega(1) = omega_avg;
    
    % Method 1: Simple integration
    for i = 2:length(theta)
        dtheta = theta(i) - theta(i-1);
        dt = dtheta / omega(i-1);  % Time step (varies with speed)
        
        % Update angular velocity
        omega(i) = omega(i-1) + alpha(i-1) * dt;
        
        % Prevent negative speeds
        if omega(i) <= 0
            warning('Negative angular velocity detected - increasing flywheel inertia may be needed');
            omega(i) = 0.1 * omega_avg;  % Set to small positive value
        end
    end
    
    % Method 2: Energy-based approach for verification
    % Change in kinetic energy = Work done by net torque
    omega_energy = zeros(size(theta));
    omega_energy(1) = omega_avg;
    
    for i = 2:length(theta)
        dtheta = theta(i) - theta(i-1);
        % Work done = torque * angular displacement
        dW = 0.5 * (T_net(i) + T_net(i-1)) * dtheta;
        
        % Change in kinetic energy: dKE = I * omega * d(omega)
        % dW = 0.5 * I * (omega_new^2 - omega_old^2)
        omega_squared = omega_energy(i-1)^2 + 2 * dW / I_flywheel;
        
        if omega_squared > 0
            omega_energy(i) = sqrt(omega_squared);
        else
            omega_energy(i) = 0.1 * omega_avg;
        end
    end
    
    % Use energy-based method as primary (more stable)
    omega = omega_energy;
    
    % Adjust to maintain correct average speed
    omega_actual_avg = mean(omega);
    omega = omega * (omega_avg / omega_actual_avg);
    
    % Convert to RPM
    rpm = omega * 60 / (2*pi);
    
    % Calculate actual coefficient of speed fluctuation
    omega_max = max(omega);
    omega_min = min(omega);
    omega_mean = mean(omega);
    Cs_actual = (omega_max - omega_min) / omega_mean;
    
    % Calculate speed characteristics
    rpm_max = max(rpm);
    rpm_min = min(rpm);
    rpm_mean = mean(rpm);
    speed_variation = (rpm_max - rpm_min);
    
    % Display dynamics results
    fprintf('Dynamic Simulation Results:\n');
    fprintf('  Speed Range: %.1f - %.1f RPM\n', rpm_min, rpm_max);
    fprintf('  Average Speed: %.1f RPM\n', rpm_mean);
    fprintf('  Speed Variation: %.1f RPM\n', speed_variation);
    fprintf('  Target Coefficient of Fluctuation: %.4f\n', params.flywheelCoefficientOfFluctuation);
    fprintf('  Actual Coefficient of Fluctuation: %.4f\n', Cs_actual);
    
    % Check if design meets requirements
    if Cs_actual > params.flywheelCoefficientOfFluctuation * 1.1
        warning('Actual coefficient of fluctuation (%.4f) exceeds target (%.4f) by >10%%', ...
                Cs_actual, params.flywheelCoefficientOfFluctuation);
        fprintf('  Consider increasing flywheel inertia\n');
    elseif Cs_actual < params.flywheelCoefficientOfFluctuation * 0.5
        fprintf('  Note: Flywheel may be oversized (Cs actual << Cs target)\n');
    else
        fprintf('  Design meets speed fluctuation requirements\n');
    end
    
    % Calculate percentage speed variation
    percent_variation = speed_variation / rpm_mean * 100;
    fprintf('  Speed Variation: %.1f%% of mean\n', percent_variation);
    fprintf('\n');
    
    % Validate results
    if any(omega <= 0)
        error('Non-positive angular velocity detected in simulation');
    end
    
    if Cs_actual > 0.2
        warning('Very high speed fluctuation (>20%%) - check flywheel sizing');
    end
end