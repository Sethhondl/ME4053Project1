%% Stirling Engine Cycle Problem

clear; clc; close all;

%% Functions

function pistonPosition = calculatePistonPosition(crankAngle, crankLength, rodLength)
    % Position is relative to bottom dead center (BDC)
    beta = asin(crankLength * sin(crankAngle) / rodLength);
    pistonPosition = rodLength * cos(beta) - crankLength * cos(crankAngle);
end

function coldVol = calculateColdVolume(crankAngle, params)
    % Calculate piston positions
    powerPistonPos = calculatePistonPosition(crankAngle, params.powerCrankLength, params.powerRodLength);
    displacerPos = calculatePistonPosition(crankAngle + params.phaseShift, params.displacerCrankLength, params.displacerRodLength);
    
    % Calculate cold side height
    % Distance between displacer and power piston, minus powerPinToPistonTop, minus half displacer height
    coldVol.height = (displacerPos - powerPistonPos) - params.powerPinToPistonTop - (params.displacerHeight / 2);
    
    % Calculate cold volume
    coldVol.volume = params.cylinderCrossSectionalArea * coldVol.height;
    
    % Ensure volume is non-negative
    coldVol.volume = max(coldVol.volume, 0);
end

function hotVol = calculateHotVolume(crankAngle, params)
%CALCULATEHOTVOLUME Calculate hot side volume
%   Inputs:
%       crankAngle - Crank angle in radians
%       params - Structure containing engine parameters
%   Output:
%       hotVol - Structure containing hot volume data
%           hotVol.volume - Hot side volume in m³
%           hotVol.height - Hot side height in m

    % Calculate piston positions
    powerPistonPos = calculatePistonPosition(crankAngle, params.powerCrankLength, params.powerRodLength);
    displacerPos = calculatePistonPosition(crankAngle + params.phaseShift, params.displacerCrankLength, params.displacerRodLength);
    
    % Calculate cold side height (same as in calculateColdVolume)
    coldHeight = (displacerPos - powerPistonPos) - params.powerPinToPistonTop - (params.displacerHeight / 2);
    
    % Calculate hot side height
    % Hot height = ColdHotHeight - cold height
    hotVol.height = params.totalCylinderHeight - .5 * params.displacerHeight - displacerPos;
    
    % Calculate hot volume
    hotVol.volume = params.cylinderCrossSectionalArea * hotVol.height;
    
    % Ensure volume is non-negative
    hotVol.volume = max(hotVol.volume, 0);
end

function schmidt = calculateSchmidtAnalysis(crankAngle, params)
%CALCULATESCHMIDTANALYSIS Calculate Schmidt analysis for Stirling engine
%   Inputs:
%       crankAngle - Crank angle in radians
%       params - Structure containing engine parameters
%   Output:
%       schmidt - Structure containing Schmidt analysis results
%           schmidt.pressure - Instantaneous pressure in Pa
%           schmidt.totalMass - Total mass in system in kg
%           schmidt.coldVolume - Cold side volume in m³
%           schmidt.hotVolume - Hot side volume in m³
%           schmidt.regeneratorVolume - Regenerator volume in m³

    % Calculate volumes at current crank angle
    coldVol = calculateColdVolume(crankAngle, params);
    hotVol = calculateHotVolume(crankAngle, params);
    
    % Extract volumes
    V_c = coldVol.volume;           % Cold side volume
    V_h = hotVol.volume;            % Hot side volume
    V_reg = params.regeneratorVolume; % Regenerator volume (constant)
    
    % Extract temperatures
    T_c = params.coldTemperature;   % Cold side temperature
    T_h = params.hotTemperature;    % Hot side temperature
    T_r = params.regeneratorTemperature; % Regenerator temperature
    
    % Extract gas constant
    R = params.gasConstant;         % Gas constant
    
    % Calculate total mass using ideal gas law at BDC condition
    % At BDC, we know the pressure and can calculate total mass
    V_comp_bdc = calculateColdVolume(0, params).volume;  % Cold volume at BDC
    V_exp_bdc = calculateHotVolume(0, params).volume;    % Hot volume at BDC
    P_bdc = params.pressureAtBDC;                        % Pressure at BDC
    
    % Apply ideal gas law at BDC condition
    denominator_bdc = V_comp_bdc/T_c + V_reg/T_r + V_exp_bdc/T_h;
    m_total = P_bdc * denominator_bdc / R;  % kg
    
    % Schmidt equation: P = m_total * R / (V_c/T_c + V_r/T_r + V_h/T_h)
    denominator = V_c/T_c + V_reg/T_r + V_h/T_h;
    P = (m_total * R) / denominator;  % Pa
    
    % Store results in structure
    schmidt.pressure = P;
    schmidt.totalMass = m_total;
    schmidt.coldVolume = V_c;
    schmidt.hotVolume = V_h;
    schmidt.regeneratorVolume = V_reg;
end

function torque = calculateTorque(crankAngle, params)
    % Pressure from Schmidt at this crank angle (make sure its absolute)
    schmidt = calculateSchmidtAnalysis(crankAngle, params);
    P      = schmidt.pressure;

    % Geometry / constants
    r  = params.powerCrankLength;
    l  = params.powerRodLength;
    A  = params.cylinderCrossSectionalArea;   % same as cylinderArea you use elsewhere
    Patm = params.atmosphericPressure;

    % Net axial force on power piston
    Fp = (P - Patm) * A;

    % Rod obliquity
    sb = (r/l) * sin(crankAngle);             % sinβ

    cb = sqrt(1 - sb.^2);                 % cosβ > 0
    torque.power = Fp * r * sin(crankAngle) ./ cb;

    torque.displacer = 0;                     % equal pressure both sides, zero rod area
    torque.total     = torque.power;
    torque.mean      = torque.total;          % single angle → same as total
end


function flywheel = sizeFlywheel(theta, params)
%SIZEFLYWHEEL Calculate required flywheel dimensions based on energy fluctuation
%   Inputs:
%       theta - Crank angle array in radians
%       params - Structure containing engine parameters
%   Output:
%       flywheel - Structure containing flywheel sizing results
%           flywheel.outerDiameter - Outer diameter in m
%           flywheel.requiredInertia - Required moment of inertia in kg·m²
%           flywheel.mass - Flywheel mass in kg
%           flywheel.energyFluctuation - Energy fluctuation in J
%           flywheel.innerDiameter - Inner diameter in m
%           flywheel.volume - Flywheel volume in m³

    % Extract parameters
    omega_avg = params.averageRPM * 2*pi/60;  % Convert RPM to rad/s
    Cs = params.flywheelCoefficientOfFluctuation;
    w = params.flywheelWidth;
    t = params.flywheelRimThickness;
    rho = params.flywheelMaterialDensity;
    
    % Calculate torque for entire cycle
    T_total = zeros(size(theta));
    for i = 1:length(theta)
        torque = calculateTorque(theta(i), params);
        T_total(i) = torque.total;
    end
    
    T_mean = mean(T_total);
    T_deviation = T_total - T_mean;
    
    % Calculate energy fluctuation
    energy_variation = zeros(size(theta));
    for i = 2:length(theta)
        dtheta = theta(i) - theta(i-1);
        energy_variation(i) = energy_variation(i-1) + ...
                              0.5 * (T_deviation(i) + T_deviation(i-1)) * dtheta;
    end
    
    E_max = max(energy_variation);
    E_min = min(energy_variation);
    energy_fluctuation = E_max - E_min;
    
    % Required moment of inertia
    I_required = energy_fluctuation / (Cs * omega_avg^2);
    
    % Calculate flywheel dimensions using iterative method
    r_outer_guess = (I_required / (2 * pi * rho * w * t))^(1/3);
    
    for iter = 1:10
        r_inner = r_outer_guess - t;
        if r_inner <= 0
            error('Flywheel thickness too large for required inertia');
        end
        
        V = pi * w * (r_outer_guess^2 - r_inner^2);
        m = rho * V;
        I_calc = 0.5 * m * (r_outer_guess^2 + r_inner^2);
        
        error_ratio = I_required / I_calc;
        r_outer_guess = r_outer_guess * error_ratio^(1/3);
        
        if abs(I_calc - I_required) / I_required < 0.001
            break;
        end
    end
    
    r_outer = r_outer_guess;
    D_outer = 2 * r_outer;
    r_inner = r_outer - t;
    D_inner = 2 * r_inner;
    
    V_flywheel = pi * w * (r_outer^2 - r_inner^2);
    mass = rho * V_flywheel;
    
    % Store results in structure
    flywheel.outerDiameter = D_outer;
    flywheel.innerDiameter = D_inner;
    flywheel.requiredInertia = I_required;
    flywheel.mass = mass;
    flywheel.energyFluctuation = energy_fluctuation;
    flywheel.volume = V_flywheel;
    
    % Check against maximum diameter limit
    if D_outer > params.maximumFlywheelDiameter
        warning('Flywheel diameter (%.2f m) exceeds maximum limit (%.2f m)', ...
                D_outer, params.maximumFlywheelDiameter);
    end
end

function dynamics = simulateDynamics(theta, I_flywheel, params)
%SIMULATEDYNAMICS Simulate angular velocity variation with flywheel
%   Inputs:
%       theta - Crank angle array in radians
%       I_flywheel - Flywheel moment of inertia in kg·m²
%       params - Structure containing engine parameters
%   Output:
%       dynamics - Structure containing dynamics simulation results
%           dynamics.angularVelocity - Angular velocity in rad/s
%           dynamics.angularAcceleration - Angular acceleration in rad/s²
%           dynamics.rpm - Angular velocity in RPM
%           dynamics.coefficientOfFluctuation - Actual coefficient of fluctuation
%           dynamics.netTorque - Net torque (engine - load) in N·m
%           dynamics.loadTorque - Load torque in N·m

    % Extract parameters
    omega_avg = params.averageRPM * 2*pi/60;  % Convert RPM to rad/s
    
    % Calculate torque for entire cycle
    T_total = zeros(size(theta));
    for i = 1:length(theta)
        torque = calculateTorque(theta(i), params);
        T_total(i) = torque.total;
    end
    
    T_mean = mean(T_total);
    T_load = T_mean;  % Load torque equals mean engine torque
    T_net = T_total - T_load;  % Net torque available for acceleration
    
    % Calculate angular acceleration
    alpha = T_net / I_flywheel;
    
    % Energy-based approach for velocity calculation
    omega = zeros(size(theta));
    omega(1) = omega_avg;  % Start at average velocity
    
    for i = 2:length(theta)
        dtheta = theta(i) - theta(i-1);
        dW = 0.5 * (T_net(i) + T_net(i-1)) * dtheta;  % Work done by net torque
        omega_squared = omega(i-1)^2 + 2 * dW / I_flywheel;
        
        if omega_squared > 0
            omega(i) = sqrt(omega_squared);
        else
            omega(i) = 0.1 * omega_avg;  % Minimum velocity to avoid zero
        end
    end
    
    % Adjust to maintain correct average speed
    omega_actual_avg = mean(omega);
    omega = omega * (omega_avg / omega_actual_avg);
    
    % Convert to RPM
    rpm = omega * 60 / (2*pi);
    
    % Calculate actual coefficient of fluctuation
    omega_max = max(omega);
    omega_min = min(omega);
    omega_mean = mean(omega);
    Cs_actual = (omega_max - omega_min) / omega_mean;
    
    % Store results in structure
    dynamics.angularVelocity = omega;
    dynamics.angularAcceleration = alpha;
    dynamics.rpm = rpm;
    dynamics.coefficientOfFluctuation = Cs_actual;
    dynamics.netTorque = T_net;
    dynamics.loadTorque = T_load;
    
    % Error checking
    if any(omega <= 0)
        error('Non-positive angular velocity detected');
    end
    
    if Cs_actual > params.flywheelCoefficientOfFluctuation * 1.1
        warning('Actual coefficient of fluctuation (%.4f) exceeds target (%.4f)', ...
                Cs_actual, params.flywheelCoefficientOfFluctuation);
    end
end



%% Prescribed Parameters

% ============== GEOMETRY PARAMETERS ==============

% Power Piston
params.powerCrankLength = 0.025;              % m - Crank radius
params.powerRodLength = 0.075;                % m - Connecting rod length
params.powerPinToPistonTop = 0.005;           % m - Pin to piston distance

% Displacer
params.displacerCrankLength = 0.020;          % m - Crank radius
params.displacerRodLength = 0.140;            % m - Connecting rod length
params.displacerVolume = 4e-5;                % m³ - Displacer volume

% Cylinder
params.cylinderBore = 0.050;                  % m - Bore diameter



% ============== OPERATING CONDITIONS ==============

% Kinematics
params.phaseShift = pi/2;                     % rad - Phase angle (90°)
params.averageRPM = 650;                      % RPM - Operating speed

% Thermodynamics
params.hotTemperature = 900;                  % K - Hot space temperature
params.coldTemperature = 300;                 % K - Cold space temperature


% Pressure
params.pressureAtBDC = 500e3;                 % Pa - Pressure at bottom dead center
params.atmosphericPressure = 101.3e3;         % Pa - Atmospheric pressure


% ============== DEAD VOLUMES ==============

params.compressionRatio = 1.7;                % Given compression ratio
params.regeneratorVolume = 2e-5;              % m³ - Given regenerator volume


% ============== WORKING FLUID (AIR) ==============

params.gasConstant = 287;                     % J/(kg·K) - Specific gas constant
params.gasGamma = 1.4;                        % - Heat capacity ratio
params.gasName = 'Air';


% ============== FLYWHEEL SPECIFICATIONS ==============

params.flywheelWidth = 0.025;                 % m - Flywheel width
params.flywheelRimThickness = 0.050;          % m - Rim thickness
params.flywheelMaterialDensity = 7750;        % kg/m³ - 304 Stainless Steel
params.flywheelCoefficientOfFluctuation = 0.003;  % - Target Cs value

% ============== SIMULATION PARAMETERS ==============

params.simulationPointsPerCycle = 360;        % Points per cycle
params.simulationCycles = 3;                  % Number of cycles
params.simulationTolerance = 1e-6;            % Convergence tolerance
params.maximumFlywheelDiameter = 2.0;         % m - Maximum allowable flywheel diameter


%% Assisting Calculations

% Cylinder cross-sectional area
params.cylinderCrossSectionalArea = pi/4*(params.cylinderBore)^2;

% Displacer GEOMETRY
params.displacerHeight = params.displacerVolume / params.cylinderCrossSectionalArea;

% Calculate power piston positions at BDC and TDC
params.powerPistonPosBDC = calculatePistonPosition(0, params.powerCrankLength, params.powerRodLength);  % 0 radians = BDC
params.powerPistonPosTDC = calculatePistonPosition(pi, params.powerCrankLength, params.powerRodLength);  % π radians = TDC

% Calculate swept volume (volume displaced by power piston)
params.powerSweptVolume = params.cylinderCrossSectionalArea * (params.powerPistonPosTDC - params.powerPistonPosBDC);

% Calculate total volume at BDC (9.23 Pg 2) (Also V_MAX)
params.totalVolumeBDC = params.regeneratorVolume - params.displacerVolume + ( params.compressionRatio * params.powerSweptVolume ) / ( params.compressionRatio - 1 ) ;

%Total cylinder height (crank pin to cylinder roof)
params.ColdHotHeight = params.totalVolumeBDC / params.cylinderCrossSectionalArea;
params.totalCylinderHeight = params.ColdHotHeight + params.displacerHeight + params.powerPinToPistonTop + params.powerRodLength -params.powerCrankLength;

% Regenerator temperature
params.regeneratorTemperature = (params.hotTemperature + params.coldTemperature) / 2;  % K












%% Test Animation - Stirling Engine Cycle Analysis

% Create crank angle array for one complete cycle
theta = linspace(0, 2*pi, params.simulationPointsPerCycle);

% Initialize arrays for storage
powerPistonPos = zeros(size(theta));
displacerPos = zeros(size(theta));
totalVolume = zeros(size(theta));
hotVolume = zeros(size(theta));
coldVolume = zeros(size(theta));
regeneratorVolume = params.regeneratorVolume * ones(size(theta));
pressure = zeros(size(theta));
totalTorque = zeros(size(theta));
powerTorque = zeros(size(theta));

% Calculate all data for each crank angle
fprintf('Calculating cycle data...\n');
for i = 1:length(theta)
    % Calculate piston positions
    powerPistonPos(i) = calculatePistonPosition(theta(i), params.powerCrankLength, params.powerRodLength);
    displacerPos(i) = calculatePistonPosition(theta(i) + params.phaseShift, params.displacerCrankLength, params.displacerRodLength);
    
    % Calculate volumes
    coldVol = calculateColdVolume(theta(i), params);
    hotVol = calculateHotVolume(theta(i), params);
    
    coldVolume(i) = coldVol.volume;
    hotVolume(i) = hotVol.volume;
    totalVolume(i) = coldVolume(i) + hotVolume(i) + regeneratorVolume(i);
    
    % Calculate Schmidt analysis
    schmidt = calculateSchmidtAnalysis(theta(i), params);
    pressure(i) = schmidt.pressure;
    
    % Calculate torque
    torque = calculateTorque(theta(i), params);
    totalTorque(i) = torque.total;
    powerTorque(i) = torque.power;
end

% Calculate flywheel sizing and dynamics
fprintf('Calculating flywheel size and dynamics...\n');
flywheel = sizeFlywheel(theta, params);
dynamics = simulateDynamics(theta, flywheel.requiredInertia, params);

% Calculate cycle averages and limits
meanPressure = mean(pressure);
maxPressure = max(pressure);
minPressure = min(pressure);
meanTorque = mean(totalTorque);
meanAngularVelocity = mean(dynamics.rpm);

% Create comprehensive analysis figure
figure('Name', 'Stirling Engine Cycle Analysis', 'Position', [50, 50, 1800, 1000]);

% Subplot 1: Piston Positions vs Crank Angle
subplot(2,3,1);
plot(theta*180/pi, powerPistonPos*1000, 'b-', 'LineWidth', 2, 'DisplayName', 'Power Piston');
hold on;
plot(theta*180/pi, displacerPos*1000, 'r-', 'LineWidth', 2, 'DisplayName', 'Displacer');
xlabel('Crank Angle (degrees)');
ylabel('Position (mm)');
title('Piston Positions vs Crank Angle');
legend('Location', 'best');
grid on;

% Subplot 2: Volume vs Crank Angle
subplot(2,3,2);
plot(theta*180/pi, totalVolume*1e6, 'k-', 'LineWidth', 2, 'DisplayName', 'Total Volume');
hold on;
plot(theta*180/pi, hotVolume*1e6, 'r-', 'LineWidth', 2, 'DisplayName', 'Hot Volume');
plot(theta*180/pi, coldVolume*1e6, 'b-', 'LineWidth', 2, 'DisplayName', 'Cold Volume');
plot(theta*180/pi, regeneratorVolume*1e6, 'g-', 'LineWidth', 2, 'DisplayName', 'Regenerator Volume');
xlabel('Crank Angle (degrees)');
ylabel('Volume (cm³)');
title('Volumes vs Crank Angle');
legend('Location', 'best');
grid on;

% Subplot 3: Pressure vs Crank Angle with constant pressure lines
subplot(2,3,3);
plot(theta*180/pi, pressure/1000, 'k-', 'LineWidth', 2, 'DisplayName', 'Instantaneous Pressure');
hold on;
plot(theta*180/pi, meanPressure/1000*ones(size(theta)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Mean Pressure');
plot(theta*180/pi, maxPressure/1000*ones(size(theta)), 'g--', 'LineWidth', 1, 'DisplayName', 'Max Pressure');
plot(theta*180/pi, minPressure/1000*ones(size(theta)), 'b--', 'LineWidth', 1, 'DisplayName', 'Min Pressure');
xlabel('Crank Angle (degrees)');
ylabel('Pressure (kPa)');
title('Pressure vs Crank Angle');
legend('Location', 'best');
grid on;

% Subplot 4: P-V Diagram with constant pressure/temperature lines
subplot(2,3,4);
plot(totalVolume*1e6, pressure/1000, 'k-', 'LineWidth', 2, 'DisplayName', 'Cycle');
hold on;
% Add constant pressure lines
plot([min(totalVolume)*1e6, max(totalVolume)*1e6], [meanPressure/1000, meanPressure/1000], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Mean Pressure');
plot([min(totalVolume)*1e6, max(totalVolume)*1e6], [maxPressure/1000, maxPressure/1000], 'g--', 'LineWidth', 1, 'DisplayName', 'Max Pressure');
plot([min(totalVolume)*1e6, max(totalVolume)*1e6], [minPressure/1000, minPressure/1000], 'b--', 'LineWidth', 1, 'DisplayName', 'Min Pressure');
xlabel('Total Volume (cm³)');
ylabel('Pressure (kPa)');
title('P-V Diagram (Schmidt Analysis)');
legend('Location', 'best');
grid on;
axis tight;

% Subplot 5: Torque vs Crank Angle with average
subplot(2,3,5);
plot(theta*180/pi, totalTorque, 'k-', 'LineWidth', 2, 'DisplayName', 'Total Torque');
hold on;
plot(theta*180/pi, powerTorque, 'b-', 'LineWidth', 2, 'DisplayName', 'Power Piston Torque');
plot(theta*180/pi, meanTorque*ones(size(theta)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Mean Torque');
xlabel('Crank Angle (degrees)');
ylabel('Torque (N·m)');
title('Torque vs Crank Angle');
legend('Location', 'best');
grid on;

% Subplot 6: Angular Velocity vs Crank Angle with average
subplot(2,3,6);
plot(theta*180/pi, dynamics.rpm, 'k-', 'LineWidth', 2, 'DisplayName', 'Angular Velocity');
hold on;
plot(theta*180/pi, meanAngularVelocity*ones(size(theta)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Mean Velocity');
plot(theta*180/pi, params.averageRPM*ones(size(theta)), 'g--', 'LineWidth', 1, 'DisplayName', 'Target Velocity');
xlabel('Crank Angle (degrees)');
ylabel('Angular Velocity (RPM)');
title('Angular Velocity vs Crank Angle');
legend('Location', 'best');
grid on;

% Display comprehensive analysis results
fprintf('\n=== STIRLING ENGINE CYCLE ANALYSIS RESULTS ===\n');
fprintf('Engine Configuration:\n');
fprintf('  Power Piston Crank Length: %.1f mm\n', params.powerCrankLength*1000);
fprintf('  Power Piston Rod Length: %.1f mm\n', params.powerRodLength*1000);
fprintf('  Displacer Crank Length: %.1f mm\n', params.displacerCrankLength*1000);
fprintf('  Displacer Rod Length: %.1f mm\n', params.displacerRodLength*1000);
fprintf('  Phase Shift: %.1f degrees\n', params.phaseShift*180/pi);
fprintf('  Cylinder Bore: %.1f mm\n', params.cylinderBore*1000);

fprintf('\nThermodynamic Conditions:\n');
fprintf('  Hot Temperature: %.0f K (%.0f°C)\n', params.hotTemperature, params.hotTemperature-273.15);
fprintf('  Cold Temperature: %.0f K (%.0f°C)\n', params.coldTemperature, params.coldTemperature-273.15);
fprintf('  Regenerator Temperature: %.0f K (%.0f°C)\n', params.regeneratorTemperature, params.regeneratorTemperature-273.15);
fprintf('  Pressure at BDC: %.1f kPa\n', params.pressureAtBDC/1000);

fprintf('\nPiston Motion Analysis:\n');
fprintf('  Power Piston Stroke: %.2f mm\n', (max(powerPistonPos) - min(powerPistonPos))*1000);
fprintf('  Displacer Stroke: %.2f mm\n', (max(displacerPos) - min(displacerPos))*1000);
fprintf('  Power Piston Max Position: %.2f mm\n', max(powerPistonPos)*1000);
fprintf('  Displacer Max Position: %.2f mm\n', max(displacerPos)*1000);

fprintf('\nVolume Analysis:\n');
fprintf('  Total Volume Range: %.2f - %.2f cm³\n', min(totalVolume)*1e6, max(totalVolume)*1e6);
fprintf('  Compression Ratio: %.2f\n', max(totalVolume)/min(totalVolume));
fprintf('  Hot Volume Range: %.2f - %.2f cm³\n', min(hotVolume)*1e6, max(hotVolume)*1e6);
fprintf('  Cold Volume Range: %.2f - %.2f cm³\n', min(coldVolume)*1e6, max(coldVolume)*1e6);
fprintf('  Regenerator Volume: %.2f cm³ (constant)\n', regeneratorVolume(1)*1e6);

fprintf('\nPressure Analysis (Schmidt):\n');
fprintf('  Pressure Range: %.2f - %.2f kPa\n', min(pressure)/1000, max(pressure)/1000);
fprintf('  Mean Pressure: %.2f kPa\n', meanPressure/1000);
fprintf('  Pressure Ratio: %.2f\n', max(pressure)/min(pressure));

fprintf('\nTorque Analysis:\n');
fprintf('  Total Torque Range: %.3f - %.3f N·m\n', min(totalTorque), max(totalTorque));
fprintf('  Mean Total Torque: %.3f N·m\n', meanTorque);
fprintf('  Power Piston Torque Range: %.3f - %.3f N·m\n', min(powerTorque), max(powerTorque));
fprintf('  Mean Power Piston Torque: %.3f N·m\n', mean(powerTorque));

fprintf('\nFlywheel Analysis:\n');
fprintf('  Required Moment of Inertia: %.4f kg·m²\n', flywheel.requiredInertia);
fprintf('  Flywheel Outer Diameter: %.3f m (%.1f mm)\n', flywheel.outerDiameter, flywheel.outerDiameter*1000);
fprintf('  Flywheel Inner Diameter: %.3f m (%.1f mm)\n', flywheel.innerDiameter, flywheel.innerDiameter*1000);
fprintf('  Flywheel Mass: %.2f kg\n', flywheel.mass);
fprintf('  Flywheel Volume: %.6f m³ (%.1f cm³)\n', flywheel.volume, flywheel.volume*1e6);
fprintf('  Energy Fluctuation: %.3f J\n', flywheel.energyFluctuation);

fprintf('\nDynamics Analysis:\n');
fprintf('  Angular Velocity Range: %.1f - %.1f RPM\n', min(dynamics.rpm), max(dynamics.rpm));
fprintf('  Mean Angular Velocity: %.1f RPM\n', meanAngularVelocity);
fprintf('  Target Angular Velocity: %.1f RPM\n', params.averageRPM);
fprintf('  Actual Coefficient of Fluctuation: %.4f\n', dynamics.coefficientOfFluctuation);
fprintf('  Target Coefficient of Fluctuation: %.4f\n', params.flywheelCoefficientOfFluctuation);
fprintf('  Net Torque Range: %.3f - %.3f N·m\n', min(dynamics.netTorque), max(dynamics.netTorque));
fprintf('  Load Torque: %.3f N·m\n', dynamics.loadTorque);

fprintf('\nPerformance Summary:\n');
fprintf('  Engine Power Output: %.2f W (estimated from mean torque)\n', meanTorque * params.averageRPM * 2*pi/60);
fprintf('  Flywheel Effectiveness: %.1f%% (target Cs: %.4f, actual Cs: %.4f)\n', ...
        (1 - dynamics.coefficientOfFluctuation/params.flywheelCoefficientOfFluctuation)*100, ...
        params.flywheelCoefficientOfFluctuation, dynamics.coefficientOfFluctuation);
fprintf('  Cycle Completeness: %.1f%% (pressure returns to within 1%% of start)\n', ...
        (1 - abs(pressure(end) - pressure(1))/pressure(1))*100);

fprintf('\n===============================================\n');