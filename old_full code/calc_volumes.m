function [totalVolume, expansionVolume, compressionVolume, powerPistonPosition, displacerPosition] = calc_volumes(crankAngle, params)
    % CALC_VOLUMES Calculate instantaneous volumes in Stirling engine
    %
    % Inputs:
    %   crankAngle - crank angle (radians)
    %   params - engine parameters structure
    %
    % Outputs:
    %   totalVolume - total gas volume (m^3)
    %   expansionVolume - hot space volume (m^3)
    %   compressionVolume - cold space volume (m^3)
    %   powerPistonPosition - position of power piston (m)
    %   displacerPosition - position of displacer (m)

    % Calculate piston positions using crank-slider kinematics
    powerPistonPos = calculatePistonPosition(crankAngle, params.powerCrankLength, params.powerRodLength);
    displacerPos = calculatePistonPosition(crankAngle + params.phaseShift, params.displacerCrankLength, params.displacerRodLength);

    % Calculate volumes based on physical geometry
    coldVol = calculateColdVolume(crankAngle, params);
    hotVol = calculateHotVolume(crankAngle, params);

    % Assign outputs
    compressionVolume = coldVol.volume;
    expansionVolume = hotVol.volume;
    totalVolume = compressionVolume + expansionVolume + params.regeneratorVolume;

    % Store piston positions
    powerPistonPosition = powerPistonPos;
    displacerPosition = displacerPos;

    % Verify volume conservation
    volumeConservationCheck = abs(totalVolume - (compressionVolume + expansionVolume + params.regeneratorVolume));
    if any(volumeConservationCheck > 1e-10)
        error('Volume conservation violated');
    end

    % Check for negative volumes
    if any(compressionVolume < 0) || any(expansionVolume < 0)
        error('Negative volume detected');
    end

    % Display volume analysis on first call
    persistent displayed;
    if isempty(displayed)
        maximumVolume = max(totalVolume);
        minimumVolume = min(totalVolume);
        actualCompressionRatio = maximumVolume / minimumVolume;

        fprintf('Volume Analysis:\n');
        fprintf('  Max Volume: %.2f cm³\n', maximumVolume * 1e6);
        fprintf('  Min Volume: %.2f cm³\n', minimumVolume * 1e6);
        fprintf('  Compression Ratio: %.2f\n', actualCompressionRatio);
        fprintf('  Swept Volume (Power): %.2f cm³\n', params.powerSweptVolume * 1e6);
        fprintf('  Dead Volume (Regenerator): %.2f cm³\n', params.regeneratorVolume * 1e6);
        fprintf('\n');

        displayed = true;
    end
end

function pistonPosition = calculatePistonPosition(crankAngle, crankLength, rodLength)
    % Calculate piston position using crank-slider kinematics
    % Position is relative to bottom dead center (BDC)
    beta = asin(crankLength * sin(crankAngle) / rodLength);
    pistonPosition = rodLength * cos(beta) - crankLength * cos(crankAngle);
end

function coldVol = calculateColdVolume(crankAngle, params)
    % Calculate cold volume based on beta engine geometry
    % Cold volume = space between displacer bottom and power piston top

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
    % Calculate hot volume based on beta engine geometry
    % Hot volume = space between displacer top and cylinder head

    % Calculate piston positions
    powerPistonPos = calculatePistonPosition(crankAngle, params.powerCrankLength, params.powerRodLength);
    displacerPos = calculatePistonPosition(crankAngle + params.phaseShift, params.displacerCrankLength, params.displacerRodLength);

    % Calculate cold side height (same as in calculateColdVolume)
    coldHeight = (displacerPos - powerPistonPos) - params.powerPinToPistonTop - (params.displacerHeight / 2);

    % Calculate hot side height
    % Hot height = ColdHotHeight - cold height
    hotVol.height = params.totalCylinderHeight - 0.5 * params.displacerHeight - displacerPos;

    % Calculate hot volume
    hotVol.volume = params.cylinderCrossSectionalArea * hotVol.height;

    % Ensure volume is non-negative
    hotVol.volume = max(hotVol.volume, 0);
end