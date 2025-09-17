function params = engine_parameters()
    % Beta-type Stirling Engine Parameters (1-10 kW range)

    % Power Piston
    params.powerCrankLength = 0.060;                % m
    params.powerRodLength = 0.240;                  % m

    % Displacer
    params.displacerCrankLength = 0.075;            % m
    params.displacerRodLength = 0.300;              % m
    params.displacerVolume = 0.002;                 % m^3
    params.displacerRodDiameter = 0.009;            % m

    % Cylinder
    params.cylinderBore = 0.180;                    % m
    params.cylinderArea = pi * (params.cylinderBore/2)^2;

    % Operating Parameters
    params.phaseShift = 90 * pi/180;                % radians
    params.compressionRatio = 3.5;

    % Temperatures
    params.hotTemperature = 600;                    % K
    params.coldTemperature = 350;                   % K
    params.regeneratorTemperature = (params.hotTemperature + params.coldTemperature) / 2;

    % Pressures
    params.pressureAtBDC = 1.5e6;                   % Pa
    params.atmosphericPressure = 101325;            % Pa

    % Dead Volumes
    params.deadVolumeHot = 0.00005;                 % m^3
    params.deadVolumeCold = 0.00005;                % m^3
    params.regeneratorVolume = 0.0001;              % m^3
    params.totalDeadVolume = params.deadVolumeHot + params.deadVolumeCold + params.regeneratorVolume;

    % Working Fluid (Air)
    params.gasConstant = 287;                       % J/(kg*K)
    params.gasGamma = 1.4;
    params.gasSpecificHeatConstantVolume = params.gasConstant / (params.gasGamma - 1);
    params.gasSpecificHeatConstantPressure = params.gasGamma * params.gasSpecificHeatConstantVolume;
    params.gasName = 'Air';

    % Flywheel
    params.flywheelWidth = 0.100;                   % m
    params.flywheelRimThickness = 0.040;            % m
    params.flywheelMaterialDensity = 7850;          % kg/m^3
    params.flywheelCoefficientOfFluctuation = 0.04;

    % Operating Speed
    params.averageRPM = 500;
    params.averageAngularVelocity = params.averageRPM * 2*pi/60;
    params.operatingFrequency = params.averageRPM / 60;

    % Simulation
    params.simulationPointsPerCycle = 360;
    params.simulationCycles = 3;
    params.simulationTolerance = 1e-6;

    % Limits
    params.minimumEfficiency = 0.0;
    params.maximumEfficiency = 1 - params.coldTemperature/params.hotTemperature;
    params.minimumPower = 0;
    params.maximumPower = 50000;
    params.maximumFlywheelDiameter = 2.0;
    params.minimumPressure = 0;

    % Calculated Parameters
    params.powerSweptVolume = params.cylinderArea * 2 * params.powerCrankLength;
    params.displacerSweptVolume = params.cylinderArea * 2 * params.displacerCrankLength;
    params.maximumVolume = params.powerSweptVolume + params.totalDeadVolume;
    params.minimumVolume = params.totalDeadVolume;

    calculatedCompressionRatio = params.maximumVolume / params.minimumVolume;
    if abs(calculatedCompressionRatio - params.compressionRatio) > 0.1
        params.compressionRatio = calculatedCompressionRatio;
    end
end