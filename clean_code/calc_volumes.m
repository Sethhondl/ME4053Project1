function [totalVolume, expansionVolume, compressionVolume, powerPistonPosition, displacerPosition] = calc_volumes(crankAngle, params)
    % Calculate instantaneous volumes for Beta-type Stirling engine

    powerCrankRadius = params.powerCrankLength;
    powerConnectingRodLength = params.powerRodLength;
    displacerCrankRadius = params.displacerCrankLength;
    displacerConnectingRodLength = params.displacerRodLength;
    phaseAngle = params.phaseShift;
    cylinderCrossSectionalArea = params.cylinderArea;

    if powerConnectingRodLength <= powerCrankRadius
        error('Power piston connecting rod must be longer than crank radius');
    end
    if displacerConnectingRodLength <= displacerCrankRadius
        error('Displacer connecting rod must be longer than crank radius');
    end

    % Power piston kinematics
    powerConnectingRodAngle = asin(powerCrankRadius * sin(crankAngle) / powerConnectingRodLength);
    powerPistonPositionFromCrank = powerCrankRadius * cos(crankAngle) + powerConnectingRodLength * cos(powerConnectingRodAngle);
    powerPistonAtTDC = powerCrankRadius + powerConnectingRodLength;
    powerPistonPosition = powerPistonAtTDC - powerPistonPositionFromCrank;

    % Displacer kinematics
    displacerCrankAngle = crankAngle - phaseAngle;
    displacerConnectingRodAngle = asin(displacerCrankRadius * sin(displacerCrankAngle) / displacerConnectingRodLength);
    displacerPositionFromCrank = displacerCrankRadius * cos(displacerCrankAngle) + displacerConnectingRodLength * cos(displacerConnectingRodAngle);
    displacerAtTDC = displacerCrankRadius + displacerConnectingRodLength;
    displacerPosition = displacerAtTDC - displacerPositionFromCrank;

    % Beta-type volume calculation
    totalVolume = params.totalDeadVolume + cylinderCrossSectionalArea * powerPistonPosition;
    workingGasVolume = totalVolume - params.totalDeadVolume;

    % Gas distribution based on displacer position
    displacerMaximumStroke = 2 * displacerCrankRadius;
    normalizedDisplacerPosition = displacerPosition / displacerMaximumStroke;
    volumeSplitFactor = 0.5 * (1 - cos(pi * normalizedDisplacerPosition));

    compressionVolume = params.deadVolumeCold + workingGasVolume .* volumeSplitFactor;
    expansionVolume = params.deadVolumeHot + workingGasVolume .* (1 - volumeSplitFactor);

    if any(compressionVolume < 0) || any(expansionVolume < 0)
        error('Negative volume detected');
    end
    if any(totalVolume < params.totalDeadVolume)
        error('Total volume less than dead volume');
    end
end