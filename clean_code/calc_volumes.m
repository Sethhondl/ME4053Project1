function [totalVolume, expansionVolume, compressionVolume, powerPistonPosition, displacerPosition] = calc_volumes(crankAngle, params)
    powerPistonPos = calculatePistonPosition(crankAngle, params.powerCrankLength, params.powerRodLength);
    displacerPos = calculatePistonPosition(crankAngle + params.phaseShift, params.displacerCrankLength, params.displacerRodLength);

    coldVol = calculateColdVolume(crankAngle, params);
    hotVol = calculateHotVolume(crankAngle, params);

    compressionVolume = coldVol.volume;
    expansionVolume = hotVol.volume;
    totalVolume = compressionVolume + expansionVolume + params.regeneratorVolume;

    powerPistonPosition = powerPistonPos;
    displacerPosition = displacerPos;

    volumeConservationCheck = abs(totalVolume - (compressionVolume + expansionVolume + params.regeneratorVolume));
    if any(volumeConservationCheck > 1e-10)
        error('Volume conservation violated');
    end

    if any(compressionVolume < 0) || any(expansionVolume < 0)
        error('Negative volume detected');
    end
end

function pistonPosition = calculatePistonPosition(crankAngle, crankLength, rodLength)
    beta = asin(crankLength * sin(crankAngle) / rodLength);
    pistonPosition = rodLength * cos(beta) - crankLength * cos(crankAngle);
end

function coldVol = calculateColdVolume(crankAngle, params)
    powerPistonPos = calculatePistonPosition(crankAngle, params.powerCrankLength, params.powerRodLength);
    displacerPos = calculatePistonPosition(crankAngle + params.phaseShift, params.displacerCrankLength, params.displacerRodLength);

    coldVol.height = (displacerPos - powerPistonPos) - params.powerPinToPistonTop - (params.displacerHeight / 2);

    coldVol.volume = params.cylinderCrossSectionalArea * coldVol.height;

    coldVol.volume = max(coldVol.volume, 0);
end

function hotVol = calculateHotVolume(crankAngle, params)
    powerPistonPos = calculatePistonPosition(crankAngle, params.powerCrankLength, params.powerRodLength);
    displacerPos = calculatePistonPosition(crankAngle + params.phaseShift, params.displacerCrankLength, params.displacerRodLength);

    coldHeight = (displacerPos - powerPistonPos) - params.powerPinToPistonTop - (params.displacerHeight / 2);

    hotVol.height = params.totalCylinderHeight - 0.5 * params.displacerHeight - displacerPos;

    hotVol.volume = params.cylinderCrossSectionalArea * hotVol.height;

    hotVol.volume = max(hotVol.volume, 0);
end