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












%% Test Animation - Piston Positions and Volume Analysis

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

% Calculate positions, volumes, and Schmidt analysis for each crank angle
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
end

% Create animation figure
figure('Name', 'Stirling Engine Volume Analysis', 'Position', [100, 100, 1400, 900]);

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

% Subplot 3: P-V Diagram (Schmidt Analysis)
subplot(2,3,3);
plot(totalVolume*1e6, pressure/1000, 'k-', 'LineWidth', 2);
xlabel('Total Volume (cm³)');
ylabel('Pressure (kPa)');
title('P-V Diagram (Schmidt Analysis)');
grid on;
axis tight;

% Subplot 4: Pressure vs Crank Angle
subplot(2,3,4);
plot(theta*180/pi, pressure/1000, 'k-', 'LineWidth', 2);
xlabel('Crank Angle (degrees)');
ylabel('Pressure (kPa)');
title('Pressure vs Crank Angle');
grid on;

% Subplot 5: Animated P-V Diagram
subplot(2,3,[5,6]);
h_pv = plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
hold on;
plot(totalVolume*1e6, pressure/1000, 'k-', 'LineWidth', 1, 'Color', [0.7 0.7 0.7]);
xlabel('Total Volume (cm³)');
ylabel('Pressure (kPa)');
title('Animated P-V Diagram');
grid on;
axis tight;

% Animation loop
for cycle = 1:params.simulationCycles
    for i = 1:length(theta)
        % Update P-V diagram point
        set(h_pv, 'XData', totalVolume(i)*1e6, 'YData', pressure(i)/1000);
        
        % Update title to show current cycle and angle
        title(sprintf('Animated P-V Diagram - Cycle %d/%d, Angle: %.1f°', ...
              cycle, params.simulationCycles, theta(i)*180/pi));
        
        % Pause for animation speed
        pause(0.01);
        
        % Break if figure is closed
        if ~ishandle(gcf)
            break;
        end
    end
    
    % Break if figure is closed
    if ~ishandle(gcf)
        break;
    end
end

% Display Schmidt analysis statistics
fprintf('\n=== Schmidt Analysis Results ===\n');
fprintf('Pressure Range: %.2f - %.2f kPa\n', min(pressure)/1000, max(pressure)/1000);
fprintf('Total Volume Range: %.2f - %.2f cm³\n', min(totalVolume)*1e6, max(totalVolume)*1e6);
fprintf('Compression Ratio: %.2f\n', max(totalVolume)/min(totalVolume));
fprintf('Hot Volume Range: %.2f - %.2f cm³\n', min(hotVolume)*1e6, max(hotVolume)*1e6);
fprintf('Cold Volume Range: %.2f - %.2f cm³\n', min(coldVolume)*1e6, max(coldVolume)*1e6);
fprintf('Regenerator Volume: %.2f cm³ (constant)\n', regeneratorVolume(1)*1e6);
fprintf('Power Piston Stroke: %.2f mm\n', (max(powerPistonPos) - min(powerPistonPos))*1000);
fprintf('Displacer Stroke: %.2f mm\n', (max(displacerPos) - min(displacerPos))*1000);
fprintf('Phase Shift: %.1f degrees\n', params.phaseShift*180/pi);
fprintf('===============================\n');