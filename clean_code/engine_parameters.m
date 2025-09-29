function params = engine_parameters()
    params.powerCrankLength = 0.025;
    params.powerRodLength = 0.075;
    params.powerPinToPistonTop = 0.005;

    params.displacerCrankLength = 0.020;
    params.displacerRodLength = 0.140;
    params.displacerVolume = 4e-5;
    params.displacerRodDiameter = 0.009;

    params.cylinderBore = 0.050;

    params.phaseShift = pi/2;
    params.compressionRatio = 1.7;

    params.hotTemperature = 900;
    params.coldTemperature = 300;

    params.pressureAtBDC = 500e3;
    params.atmosphericPressure = 101.3e3;

    params.regeneratorVolume = 2e-5;

    params.gasConstant = 287;
    params.gasGamma = 1.4;
    params.gasName = 'Air';

    params.flywheelWidth = 0.025;
    params.flywheelRimThickness = 0.050;
    params.flywheelMaterialDensity = 7750;
    params.flywheelCoefficientOfFluctuation = 0.003;

    params.averageRPM = 650;

    params.simulationPointsPerCycle = 360;
    params.simulationCycles = 3;
    params.simulationTolerance = 1e-6;

    params.minimumEfficiency = 0.0;
    params.minimumPower = 0;
    params.maximumPower = 50000;
    params.maximumFlywheelDiameter = 2.0;
    params.minimumPressure = 0;

    params.cylinderCrossSectionalArea = pi/4*(params.cylinderBore)^2;

    params.displacerHeight = params.displacerVolume / params.cylinderCrossSectionalArea;

    powerPistonPosBDC = params.powerRodLength * (1 - cos(asin(params.powerCrankLength * sin(0) / params.powerRodLength))) - params.powerCrankLength * cos(0) + params.powerRodLength + params.powerCrankLength;
    powerPistonPosTDC = params.powerRodLength * (1 - cos(asin(params.powerCrankLength * sin(pi) / params.powerRodLength))) - params.powerCrankLength * cos(pi) + params.powerRodLength + params.powerCrankLength;

    params.powerSweptVolume = params.cylinderCrossSectionalArea * (powerPistonPosTDC - powerPistonPosBDC);

    params.totalVolumeBDC = params.regeneratorVolume - params.displacerVolume + (params.compressionRatio * params.powerSweptVolume) / (params.compressionRatio - 1);

    params.ColdHotHeight = params.totalVolumeBDC / params.cylinderCrossSectionalArea;
    params.totalCylinderHeight = params.ColdHotHeight + params.displacerHeight + params.powerPinToPistonTop + params.powerRodLength - params.powerCrankLength;

    params.regeneratorTemperature = (params.hotTemperature + params.coldTemperature) / 2;
end