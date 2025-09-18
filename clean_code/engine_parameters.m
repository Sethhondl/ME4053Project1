function params = engine_parameters()
%ENGINE_PARAMETERS Define all Stirling engine configuration parameters
%   Returns a structure containing all engine design parameters

    %% ============== GEOMETRY PARAMETERS ==============

    % Power Piston
    params.powerCrankLength = 0.025;              % m - Crank radius
    params.powerRodLength = 0.075;                % m - Connecting rod
    params.powerPinToPistonTop = 0.005;           % m - Pin to piston distance

    % Displacer
    params.displacerCrankLength = 0.020;          % m - Crank radius
    params.displacerRodLength = 0.140;            % m - Connecting rod
    params.displacerVolume = 4e-5;                % m³ - Displacer volume
    params.displacerRodDiameter = 0.009;          % m - Rod diameter

    % Cylinder
    params.cylinderBore = 0.050;                  % m - Bore diameter
    params.cylinderArea = pi * (params.cylinderBore/2)^2;  % m² - Cross-sectional area

    % Calculated Swept Volumes
    params.powerSweptVolume = params.cylinderArea * 2 * params.powerCrankLength;      % m³
    params.displacerSweptVolume = params.cylinderArea * 2 * params.displacerCrankLength;  % m³


    %% ============== OPERATING CONDITIONS ==============

    % Kinematics
    params.phaseShift = pi/2;                     % rad - Phase angle (90°)
    params.averageRPM = 650;                      % RPM - Operating speed
    params.averageAngularVelocity = params.averageRPM * 2*pi/60;  % rad/s
    params.operatingFrequency = params.averageRPM / 60;            % Hz

    % Thermodynamics
    params.hotTemperature = 900;                  % K - Hot space temperature
    params.coldTemperature = 300;                 % K - Cold space temperature
    params.regeneratorTemperature = (params.hotTemperature + params.coldTemperature) / 2;  % K

    % Pressure
    params.pressureAtBDC = 500e3;                 % Pa - Pressure at bottom dead center
    params.atmosphericPressure = 101.3e3;         % Pa - Atmospheric pressure


    %% ============== DEAD VOLUMES (CR = 1.7) ==============

    targetCompressionRatio = 1.7;
    requiredTotalDeadVolume = params.powerSweptVolume / (targetCompressionRatio - 1);

    params.regeneratorVolume = 2e-5;              % m³ - Given regenerator volume
    remainingDeadVolume = requiredTotalDeadVolume - params.regeneratorVolume;

    params.deadVolumeHot = 0.4 * remainingDeadVolume;   % m³ - Hot space dead volume
    params.deadVolumeCold = 0.6 * remainingDeadVolume;  % m³ - Cold space dead volume
    params.totalDeadVolume = params.deadVolumeHot + params.deadVolumeCold + params.regeneratorVolume;


    %% ============== WORKING FLUID (AIR) ==============

    params.gasConstant = 287;                     % J/(kg·K) - Specific gas constant
    params.gasGamma = 1.4;                        % - Heat capacity ratio
    params.gasSpecificHeatConstantVolume = params.gasConstant / (params.gasGamma - 1);
    params.gasSpecificHeatConstantPressure = params.gasGamma * params.gasSpecificHeatConstantVolume;
    params.gasName = 'Air';


    %% ============== FLYWHEEL SPECIFICATIONS ==============

    params.flywheelWidth = 0.025;                 % m - Flywheel width
    params.flywheelRimThickness = 0.050;          % m - Rim thickness
    params.flywheelMaterialDensity = 8000;        % kg/m³ - 304 Stainless Steel
    params.flywheelCoefficientOfFluctuation = 0.003;  % - Target Cs value


    %% ============== SIMULATION PARAMETERS ==============

    params.simulationPointsPerCycle = 360;        % Points per cycle
    params.simulationCycles = 3;                  % Number of cycles
    params.simulationTolerance = 1e-6;            % Convergence tolerance


    %% ============== VALIDATION LIMITS ==============

    params.minimumEfficiency = 0.0;
    params.maximumEfficiency = 1 - params.coldTemperature/params.hotTemperature;  % Carnot limit
    params.minimumPower = 0;                      % W
    params.maximumPower = 50000;                  % W
    params.maximumFlywheelDiameter = 2.0;         % m
    params.minimumPressure = 0;                   % Pa


    %% ============== DERIVED PARAMETERS ==============

    params.compressionRatio = 1.7;                % Given compression ratio
    params.maximumVolume = params.powerSweptVolume + params.totalDeadVolume;
    params.minimumVolume = params.totalDeadVolume;

    % Verify compression ratio
    calculatedCompressionRatio = params.maximumVolume / params.minimumVolume;
    if abs(calculatedCompressionRatio - params.compressionRatio) > 0.01
        params.compressionRatio = calculatedCompressionRatio;
    end

end