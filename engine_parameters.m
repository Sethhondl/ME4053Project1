function params = engine_parameters()
    % STIRLING ENGINE PARAMETERS
    % Configurable parameters for Beta-type Stirling Engine Analysis
    % Large-scale power generation system (1-10 kW range)
    % Working Fluid: Air
    % Flywheel Material: Steel

    %% Power Piston Parameters
    params.powerCrankLength = 0.025;        % m (25 mm) - from given parameters
    params.powerRodLength = 0.075;          % m (75 mm) - from given parameters
    params.powerPinToPistonTop = 0.005;     % m (5 mm) - distance from pin to piston top

    %% Displacer Parameters
    params.displacerCrankLength = 0.020;         % m (20 mm) - from given parameters
    params.displacerRodLength = 0.140;           % m (140 mm) - from given parameters
    params.displacerVolume = 4e-5;               % m^3 (40 cm^3) - from given parameters
    params.displacerRodDiameter = 0.009;         % m (9 mm) - 5% of bore diameter for Beta-type

    %% Cylinder Parameters
    params.cylinderBore = 0.050;             % m (50 mm diameter) - from given parameters
    params.cylinderArea = pi * (params.cylinderBore/2)^2;  % m^2

    %% Operating Parameters
    params.phaseShift = pi/2;                % radians (90 degrees) - from given parameters
    params.compressionRatio = 1.7;           % dimensionless - from given parameters

    %% Temperature Conditions
    params.hotTemperature = 900;                       % K (hot temperature) - from given parameters
    params.coldTemperature = 300;                      % K (cold temperature) - from given parameters

    %% Pressure Conditions
    params.pressureAtBDC = 500e3;                     % Pa (500 kPa at bottom dead center) - from given parameters
    params.atmosphericPressure = 101.3e3;             % Pa (101.3 kPa) - from given parameters

    %% Regenerator Volume
    params.regeneratorVolume = 2e-5;             % m^3 (20 cm^3) - from given parameters

    %% Working Fluid Properties (Air)
    params.gasConstant = 287;                       % J/(kg*K) - specific gas constant for air
    params.gasGamma = 1.4;                          % heat capacity ratio
    params.gasSpecificHeatConstantVolume = params.gasConstant / (params.gasGamma - 1);  % J/(kg*K)
    params.gasSpecificHeatConstantPressure = params.gasGamma * params.gasSpecificHeatConstantVolume;       % J/(kg*K)
    params.gasName = 'Air';

    %% Flywheel Parameters
    params.flywheelWidth = 0.025;                   % m (25 mm) - from given parameters
    params.flywheelRimThickness = 0.050;            % m (50 mm rim thickness) - from given parameters
    params.flywheelMaterialDensity = 7750;          % kg/m^3 (304 stainless steel) - from given parameters
    params.flywheelCoefficientOfFluctuation = 0.003; % dimensionless (0.3%) - from given parameters

    %% Operating Speed
    params.averageRPM = 650;                     % RPM (average rotational speed) - from given parameters
    params.averageAngularVelocity = params.averageRPM * 2*pi/60;  % rad/s
    params.operatingFrequency = params.averageRPM / 60;   % Hz

    %% Simulation Parameters
    params.simulationPointsPerCycle = 360;                % number of points per cycle
    params.simulationCycles = 3;                          % number of cycles to simulate
    params.simulationTolerance = 1e-6;                    % convergence tolerance

    %% Validation Limits (for error checking)
    params.minimumEfficiency = 0.0;                       % minimum allowable efficiency
    params.maximumEfficiency = 1 - params.coldTemperature/params.hotTemperature;  % Carnot limit
    params.minimumPower = 0;                              % W
    params.maximumPower = 50000;                          % W (50 kW max)
    params.maximumFlywheelDiameter = 2.0;                 % m
    params.minimumPressure = 0;                           % Pa (must be positive)

    %% ============== DERIVED GEOMETRIC CALCULATIONS ==============
    % Cylinder cross-sectional area
    params.cylinderCrossSectionalArea = pi/4*(params.cylinderBore)^2;

    % Displacer height from volume
    params.displacerHeight = params.displacerVolume / params.cylinderCrossSectionalArea;

    % Calculate power piston positions at BDC and TDC
    powerPistonPosBDC = params.powerRodLength * (1 - cos(asin(params.powerCrankLength * sin(0) / params.powerRodLength))) - params.powerCrankLength * cos(0) + params.powerRodLength + params.powerCrankLength;
    powerPistonPosTDC = params.powerRodLength * (1 - cos(asin(params.powerCrankLength * sin(pi) / params.powerRodLength))) - params.powerCrankLength * cos(pi) + params.powerRodLength + params.powerCrankLength;

    % Calculate swept volume
    params.powerSweptVolume = params.cylinderCrossSectionalArea * (powerPistonPosTDC - powerPistonPosBDC);

    % Total volume at BDC from compression ratio
    params.totalVolumeBDC = params.regeneratorVolume - params.displacerVolume + (params.compressionRatio * params.powerSweptVolume) / (params.compressionRatio - 1);

    % Cylinder heights
    params.ColdHotHeight = params.totalVolumeBDC / params.cylinderCrossSectionalArea;
    params.totalCylinderHeight = params.ColdHotHeight + params.displacerHeight + params.powerPinToPistonTop + params.powerRodLength - params.powerCrankLength;

    % Regenerator temperature
    params.regeneratorTemperature = (params.hotTemperature + params.coldTemperature) / 2;

    %% Display Configuration Summary
    fprintf('\n=== STIRLING ENGINE CONFIGURATION ===\n');
    fprintf('Engine Type: Beta-type Stirling\n');
    fprintf('Working Fluid: %s\n', params.gasName);
    fprintf('Cylinder Bore: %.0f mm\n', params.cylinderBore * 1000);
    fprintf('Power Stroke: %.0f mm\n', params.powerCrankLength * 2 * 1000);
    fprintf('Displacer Stroke: %.0f mm\n', params.displacerCrankLength * 2 * 1000);
    fprintf('Phase Shift: %.0f degrees\n', params.phaseShift * 180/pi);
    fprintf('Compression Ratio: %.2f\n', params.compressionRatio);
    fprintf('Operating Speed: %.0f RPM\n', params.averageRPM);
    fprintf('Temperature Range: %.0f K to %.0f K\n', params.coldTemperature, params.hotTemperature);
    fprintf('Carnot Efficiency Limit: %.1f%%\n', params.maximumEfficiency * 100);
    fprintf('=====================================\n\n');
end