function [totalVolume, expansionVolume, compressionVolume, powerPistonPosition, displacerPosition] = calc_volumes(crankAngle, params)
%CALC_VOLUMES Calculate instantaneous volumes for Beta-type Stirling engine
%   Uses crank-slider kinematics to determine piston positions and volumes

    %% ========== EXTRACT PARAMETERS ==========

    powerCrankRadius = params.powerCrankLength;
    powerRodLength = params.powerRodLength;
    displacerCrankRadius = params.displacerCrankLength;
    displacerRodLength = params.displacerRodLength;
    phaseAngle = params.phaseShift;
    cylinderArea = params.cylinderArea;


    %% ========== VALIDATE GEOMETRY ==========

    if powerRodLength <= powerCrankRadius
        error('Power piston connecting rod must be longer than crank radius');
    end
    if displacerRodLength <= displacerCrankRadius
        error('Displacer connecting rod must be longer than crank radius');
    end


    %% ========== POWER PISTON KINEMATICS ==========

    % Crank-slider mechanism for power piston
    powerRodAngle = asin(powerCrankRadius * sin(crankAngle) / powerRodLength);
    powerPistonFromCrank = powerCrankRadius * cos(crankAngle) + ...
                          powerRodLength * cos(powerRodAngle);
    powerPistonAtTDC = powerCrankRadius + powerRodLength;
    powerPistonPosition = powerPistonAtTDC - powerPistonFromCrank;  % Distance from TDC


    %% ========== DISPLACER KINEMATICS ==========

    % Displacer leads by phase angle
    displacerCrankAngle = crankAngle - phaseAngle;

    % Crank-slider mechanism for displacer
    displacerRodAngle = asin(displacerCrankRadius * sin(displacerCrankAngle) / displacerRodLength);
    displacerFromCrank = displacerCrankRadius * cos(displacerCrankAngle) + ...
                        displacerRodLength * cos(displacerRodAngle);
    displacerAtTDC = displacerCrankRadius + displacerRodLength;
    displacerPosition = displacerAtTDC - displacerFromCrank;  % Distance from TDC


    %% ========== VOLUME CALCULATIONS ==========

    % Total gas volume (depends only on power piston)
    totalVolume = params.totalDeadVolume + cylinderArea * powerPistonPosition;

    % Working gas volume (excluding dead volumes)
    workingGasVolume = totalVolume - params.totalDeadVolume;

    % Gas distribution based on displacer position
    displacerStroke = 2 * displacerCrankRadius;
    normalizedPosition = displacerPosition / displacerStroke;  % 0 to 1

    % Smooth sinusoidal gas transfer
    volumeSplitFactor = 0.5 * (1 - cos(pi * normalizedPosition));

    % Calculate space volumes
    compressionVolume = params.deadVolumeCold + workingGasVolume .* volumeSplitFactor;
    expansionVolume = params.deadVolumeHot + workingGasVolume .* (1 - volumeSplitFactor);


    %% ========== VALIDATION ==========

    if any(compressionVolume < 0) || any(expansionVolume < 0)
        error('Negative volume detected - check geometry');
    end
    if any(totalVolume < params.totalDeadVolume)
        error('Total volume less than dead volume - check calculations');
    end

end