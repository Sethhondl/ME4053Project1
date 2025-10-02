function [totalVolume, expansionVolume, compressionVolume, powerPistonPosition, displacerPosition] = calc_volumes(crankAngle, params)
% CALC_VOLUMES - Calculate instantaneous volumes in beta-type Stirling engine
%
% Syntax:
%   [totalVolume, expansionVolume, compressionVolume, powerPistonPosition, displacerPosition] = calc_volumes(crankAngle, params)
%
% Description:
%   Calculates the instantaneous volumes of hot space, cold space, and total
%   gas volume based on crank angle and engine geometry. Uses proper beta
%   engine physical model where displacer shuttles gas between spaces.
%
% Inputs:
%   crankAngle - Crank angle(s) in radians (can be array)
%   params     - Engine parameters structure from engine_parameters()
%
% Outputs:
%   totalVolume          - Total gas volume (m^3)
%   expansionVolume      - Hot space volume (m^3)
%   compressionVolume    - Cold space volume (m^3)
%   powerPistonPosition  - Power piston position from BDC (m)
%   displacerPosition    - Displacer position from BDC (m)
%
% Example:
%   params = engine_parameters();
%   theta = linspace(0, 2*pi, 360);
%   [V_total, V_hot, V_cold, ~, ~] = calc_volumes(theta, params);
%
% See also: ENGINE_PARAMETERS, SCHMIDT_ANALYSIS
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