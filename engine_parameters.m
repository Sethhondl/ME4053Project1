function params = engine_parameters()
    % STIRLING ENGINE PARAMETERS
    % Configurable parameters for Beta-type Stirling Engine Analysis
    % Large-scale power generation system (1-10 kW range)
    % Working Fluid: Air
    % Flywheel Material: Steel

    %% Power Piston Parameters
    params.powerCrankLength = 0.060;        % m (60 mm) - optimized for 1-10kW
    params.powerRodLength = 0.240;          % m (240 mm)

    %% Displacer Parameters
    params.displacerCrankLength = 0.075;         % m (75 mm) - 1.25x power stroke
    params.displacerRodLength = 0.300;           % m (300 mm)
    params.displacerVolume = 0.002;              % m^3 (2000 cm^3)
    params.displacerRodDiameter = 0.009;         % m (9 mm) - 5% of bore diameter for Beta-type

    %% Cylinder Parameters
    params.cylinderBore = 0.180;             % m (180 mm diameter) - sized for 1-10kW power
    params.cylinderArea = pi * (params.cylinderBore/2)^2;  % m^2

    %% Operating Parameters
    params.phaseShift = 90 * pi/180;         % radians (90 degrees)
    params.compressionRatio = 3.5;           % dimensionless

    %% Temperature Conditions
    params.hotTemperature = 600;                       % K (hot temperature) - reduced for smaller engine
    params.coldTemperature = 350;                      % K (cold temperature)
    params.regeneratorTemperature = (params.hotTemperature + params.coldTemperature) / 2;  % K (average)

    %% Pressure Conditions
    params.pressureAtBDC = 1.5e6;                     % Pa (1.5 MPa at bottom dead center) - for target power
    params.atmosphericPressure = 101325;              % Pa (atmospheric pressure)

    %% Dead Volumes
    params.deadVolumeHot = 0.00005;              % m^3 (50 cm^3 hot space dead volume)
    params.deadVolumeCold = 0.00005;             % m^3 (50 cm^3 cold space dead volume)
    params.regeneratorVolume = 0.0001;           % m^3 (100 cm^3 regenerator dead volume)
    params.totalDeadVolume = params.deadVolumeHot + params.deadVolumeCold + params.regeneratorVolume;

    %% Working Fluid Properties (Air)
    params.gasConstant = 287;                       % J/(kg*K) - specific gas constant for air
    params.gasGamma = 1.4;                          % heat capacity ratio
    params.gasSpecificHeatConstantVolume = params.gasConstant / (params.gasGamma - 1);  % J/(kg*K)
    params.gasSpecificHeatConstantPressure = params.gasGamma * params.gasSpecificHeatConstantVolume;       % J/(kg*K)
    params.gasName = 'Air';

    %% Flywheel Parameters
    params.flywheelWidth = 0.100;                   % m (100 mm)
    params.flywheelRimThickness = 0.040;            % m (40 mm rim thickness)
    params.flywheelMaterialDensity = 7850;          % kg/m^3 (steel)
    params.flywheelCoefficientOfFluctuation = 0.04; % dimensionless (4%)

    %% Operating Speed
    params.averageRPM = 500;                     % RPM (average rotational speed)
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
    % Swept volumes
    params.powerSweptVolume = params.cylinderArea * 2 * params.powerCrankLength;  % m^3
    params.displacerSweptVolume = params.cylinderArea * 2 * params.displacerCrankLength;    % m^3

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