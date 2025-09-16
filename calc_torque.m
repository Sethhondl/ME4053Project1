function [T_total, T_power, T_disp, T_mean] = calc_torque(P, theta, x_power, x_disp, params)
    % CALC_TORQUE Calculate instantaneous torque from pressure forces
    %
    % Inputs:
    %   P - pressure array (Pa)
    %   theta - crank angle array (radians)
    %   x_power - power piston position array (m)
    %   x_disp - displacer position array (m)
    %   params - engine parameters structure
    %
    % Outputs:
    %   T_total - total torque from both pistons (N*m)
    %   T_power - torque from power piston (N*m)
    %   T_disp - torque from displacer (N*m)
    %   T_mean - mean torque over cycle (N*m)
    
    % Extract parameters
    r_p = params.power.crank_length;      % Power piston crank radius
    l_p = params.power.rod_length;        % Power piston connecting rod
    r_d = params.disp.crank_length;       % Displacer crank radius
    l_d = params.disp.rod_length;         % Displacer connecting rod
    phase = params.phase_shift;           % Phase angle
    A = params.cylinder.area;             % Piston area
    P_atm = params.P_atm;                 % Atmospheric pressure
    
    % Initialize torque arrays
    T_power = zeros(size(theta));
    T_disp = zeros(size(theta));
    
    % Calculate torque for each crank angle
    for i = 1:length(theta)
        %% Power Piston Torque
        % Force on power piston (pressure difference across piston)
        F_power = (P(i) - P_atm) * A;  % Force = pressure difference * area
        
        % Crank-slider kinematics for torque
        % Torque = Force * moment arm
        % Moment arm = r * sin(theta + beta) where beta is connecting rod angle
        
        % Power piston connecting rod angle
        sin_beta_p = r_p * sin(theta(i)) / l_p;
        cos_beta_p = sqrt(1 - sin_beta_p^2);
        
        % Torque from power piston
        % T = F * r * sin(theta + beta) / cos(beta)
        % Simplified: T = F * r * (sin(theta) + sin(theta)*cos(theta)*tan(beta))
        
        % More accurate formula considering connecting rod angle
        if abs(sin_beta_p) < 1  % Valid geometry
            T_power(i) = F_power * r_p * (sin(theta(i)) + ...
                         sin_beta_p * cos(theta(i)) / cos_beta_p);
        else
            T_power(i) = 0;  % Invalid position
        end
        
        %% Displacer Torque
        % The displacer moves gas between hot and cold spaces
        % It experiences pressure on both sides, but with different areas
        
        % Displacer crank angle (with phase shift)
        theta_disp = theta(i) - phase;
        
        % For a beta-type engine, the displacer doesn't do work directly
        % But it affects the pressure distribution
        % We'll model it as having a small pressure difference
        
        % Approximate force on displacer (pressure gradient effect)
        F_disp = P(i) * A * 0.1;  % Simplified - 10% of pressure force
        
        % Displacer connecting rod angle
        sin_beta_d = r_d * sin(theta_disp) / l_d;
        if abs(sin_beta_d) < 1
            cos_beta_d = sqrt(1 - sin_beta_d^2);
            T_disp(i) = F_disp * r_d * (sin(theta_disp) + ...
                        sin_beta_d * cos(theta_disp) / cos_beta_d);
        else
            T_disp(i) = 0;
        end
    end
    
    % Total torque is the sum of both contributions
    T_total = T_power + T_disp;
    
    % Calculate mean torque
    T_mean = mean(T_total);
    
    % Calculate torque characteristics
    T_max = max(T_total);
    T_min = min(T_total);
    T_rms = sqrt(mean(T_total.^2));  % RMS torque
    
    % Torque fluctuation
    if T_mean ~= 0
        torque_fluctuation = (T_max - T_min) / abs(T_mean);
    else
        torque_fluctuation = inf;
    end
    
    % Display torque analysis (only on first call)
    persistent displayed;
    if isempty(displayed)
        fprintf('Torque Analysis:\n');
        fprintf('  Maximum Torque: %.2f N*m\n', T_max);
        fprintf('  Minimum Torque: %.2f N*m\n', T_min);
        fprintf('  Mean Torque: %.2f N*m\n', T_mean);
        fprintf('  RMS Torque: %.2f N*m\n', T_rms);
        fprintf('  Torque Fluctuation: %.2f\n', torque_fluctuation);
        
        % Check if torque is reasonable
        expected_torque = params.limits.max_power / params.omega_avg;
        if abs(T_mean) > 2 * expected_torque
            warning('Mean torque seems high compared to expected power output');
        end
        
        % Find zero-crossing points (torque reversal)
        zero_crossings = find(diff(sign(T_total)) ~= 0);
        fprintf('  Torque Reversals per Cycle: %d\n', length(zero_crossings));
        fprintf('\n');
        
        displayed = true;
    end
    
    % Validate torque
    if all(T_total == 0)
        error('Zero torque throughout cycle - check calculations');
    end
    
    if T_mean <= 0
        warning('Non-positive mean torque - engine not producing net work');
    end
end