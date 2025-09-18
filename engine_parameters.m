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

    %% Calculate Swept Volumes (needed for dead volume calculation)
    params.powerSweptVolume = params.cylinderArea * 2 * params.powerCrankLength;  % m^3
    params.displacerSweptVolume = params.cylinderArea * 2 * params.displacerCrankLength;  % m^3

    %% Operating Parameters
    params.phaseShift = pi/2;                % radians (90 degrees) - from given parameters
    params.compressionRatio = 1.7;           % dimensionless - from given parameters

    %% Temperature Conditions
    params.hotTemperature = 900;                       % K (hot temperature) - from given parameters
    params.coldTemperature = 300;                      % K (cold temperature) - from given parameters
    params.regeneratorTemperature = (params.hotTemperature + params.coldTemperature) / 2;  % K (average)

    %% Pressure Conditions
    params.pressureAtBDC = 500e3;                     % Pa (500 kPa at bottom dead center) - from given parameters
    params.atmosphericPressure = 101.3e3;             % Pa (101.3 kPa) - from given parameters

    %% Dead Volumes (calculated to achieve compression ratio of 1.7)
    % Given compression ratio = 1.7
    % CR = V_max / V_min = (swept + dead) / dead
    % Therefore: dead = swept / (CR - 1)
    targetCompressionRatio = 1.7;  % from given parameters

    % Calculate required total dead volume
    requiredTotalDeadVolume = params.powerSweptVolume / (targetCompressionRatio - 1);

    % Given regenerator volume
    params.regeneratorVolume = 2e-5;             % m^3 (20 cm^3) - from given parameters

    % Remaining dead volume to distribute between hot and cold spaces
    remainingDeadVolume = requiredTotalDeadVolume - params.regeneratorVolume;

    % Distribute dead volume (40% hot, 60% cold typical for beta-type)
    params.deadVolumeHot = 0.4 * remainingDeadVolume;     % m^3 (hot space dead volume)
    params.deadVolumeCold = 0.6 * remainingDeadVolume;    % m^3 (cold space dead volume)
    params.totalDeadVolume = params.deadVolumeHot + params.deadVolumeCold + params.regeneratorVolume;

    %% Working Fluid Properties (Air)
    params.gasConstant = 287;                       % J/(kg*K) - specific gas constant for air
    params.gasGamma = 1.4;                          % heat capacity ratio
    params.gasSpecificHeatConstantVolume = params.gasConstant / (params.gasGamma - 1);  % J/(kg*K)
    params.gasSpecificHeatConstantPressure = params.gasGamma * params.gasSpecificHeatConstantVolume;       % J/(kg*K)
    params.gasName = 'Air';

    %% Flywheel Parameters
    params.flywheelWidth = 0.025;                   % m (25 mm) - from given parameters
    params.flywheelRimThickness = 0.050;            % m (50 mm rim thickness) - from given parameters
    params.flywheelMaterialDensity = 8000;          % kg/m^3 (304 stainless steel) - from given parameters
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

    %% Calculated Derived Parameters
    % Swept volumes already calculated above for dead volume calculation

    % For beta-type engine, compression ratio is based on power piston swept volume
    % Maximum volume when power piston at BDC
    params.maximumVolume = params.powerSweptVolume + params.totalDeadVolume;
    % Minimum volume when power piston at TDC
    params.minimumVolume = params.totalDeadVolume;

    % Check compression ratio
    calculatedCompressionRatio = params.maximumVolume / params.minimumVolume;
    if abs(calculatedCompressionRatio - params.compressionRatio) > 0.1
        % Use the calculated value
        params.compressionRatio = calculatedCompressionRatio;
    end

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