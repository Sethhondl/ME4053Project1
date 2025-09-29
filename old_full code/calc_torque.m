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
    r_p = params.powerCrankLength;      % Power piston crank radius
    l_p = params.powerRodLength;        % Power piston connecting rod
    r_d = params.displacerCrankLength;       % Displacer crank radius
    l_d = params.displacerRodLength;         % Displacer connecting rod
    phase = params.phaseShift;           % Phase angle
    A = params.cylinderArea;             % Piston area
    P_atm = params.atmosphericPressure;                 % Atmospheric pressure
    
    % Initialize torque arrays
    T_power = zeros(size(theta));
    T_disp = zeros(size(theta));
    
    % Calculate torque for each crank angle
    for i = 1:length(theta)
        %% Power Piston Torque
        % Force on power piston (pressure pushes piston down during expansion)
        % Positive force when pressure > atmospheric (expansion work)
        F_power = (P(i) - P_atm) * A;  % Force = pressure difference * area

        % Crank-slider kinematics for torque
        % Torque = Force * moment arm

        % Power piston connecting rod angle
        sin_beta_p = r_p * sin(theta(i)) / l_p;
        if abs(sin_beta_p) < 1  % Valid geometry
            cos_beta_p = sqrt(1 - sin_beta_p^2);

            % Torque from power piston using exact crank-slider formula
            % T = F * r * sin(theta) / sqrt(1 - (r/l*sin(theta))^2)
            % Note: For power stroke (0 to pi), sin(theta) > 0 and F > 0 should give T > 0
            % The negative sign corrects the convention so positive pressure gives positive torque
            T_power(i) = -F_power * r_p * sin(theta(i)) / cos_beta_p;
        else
            T_power(i) = 0;  % Invalid position
        end
        
        %% Displacer Torque
        % For beta-type engine with rod diameter assumed zero per specifications
        % The displacer experiences no net force (equal pressure on both sides)
        % Schmidt analysis assumes uniform pressure throughout engine
        % Therefore, no torque contribution from displacer

        T_disp(i) = 0;  % Zero torque as rod diameter is assumed negligible
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
        expected_torque = params.maximumPower / params.averageAngularVelocity;
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