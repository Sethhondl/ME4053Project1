function [totalVolume, expansionVolume, compressionVolume, powerPistonPosition, displacerPosition] = calc_volumes(crankAngle, params)
    % CALC_VOLUMES Calculate instantaneous volumes for Beta-type Stirling engine
    %
    % Inputs:
    %   crankAngle - crank angle array (radians)
    %   params - engine parameters structure
    %
    % Outputs:
    %   totalVolume - total gas volume (m^3)
    %   expansionVolume - expansion space volume (m^3)
    %   compressionVolume - compression space volume (m^3)
    %   powerPistonPosition - power piston position (m)
    %   displacerPosition - displacer position (m)

    % Extract parameters with descriptive names
    powerCrankRadius = params.powerCrankLength;
    powerConnectingRodLength = params.powerRodLength;
    displacerCrankRadius = params.displacerCrankLength;
    displacerConnectingRodLength = params.displacerRodLength;
    phaseAngle = params.phaseShift;
    cylinderCrossSectionalArea = params.cylinderArea;

    % Validate crank-slider geometry
    if powerConnectingRodLength <= powerCrankRadius
        error('Power piston connecting rod must be longer than crank radius');
    end
    if displacerConnectingRodLength <= displacerCrankRadius
        error('Displacer connecting rod must be longer than crank radius');
    end

    % Calculate piston positions using crank-slider kinematics
    % Position measured from TDC (Top Dead Center)

    % Power piston position (compression space)
    powerConnectingRodAngle = asin(powerCrankRadius * sin(crankAngle) / powerConnectingRodLength);
    powerPistonPositionFromCrank = powerCrankRadius * cos(crankAngle) + powerConnectingRodLength * cos(powerConnectingRodAngle);
    powerPistonAtTDC = powerCrankRadius + powerConnectingRodLength;  % Position at TDC
    powerPistonPosition = powerPistonAtTDC - powerPistonPositionFromCrank;  % Distance from TDC

    % Displacer position (with phase shift)
    displacerCrankAngle = crankAngle - phaseAngle;  % Displacer leads by phase angle
    displacerConnectingRodAngle = asin(displacerCrankRadius * sin(displacerCrankAngle) / displacerConnectingRodLength);
    displacerPositionFromCrank = displacerCrankRadius * cos(displacerCrankAngle) + displacerConnectingRodLength * cos(displacerConnectingRodAngle);
    displacerAtTDC = displacerCrankRadius + displacerConnectingRodLength;
    displacerPosition = displacerAtTDC - displacerPositionFromCrank;

    % Calculate instantaneous volumes for beta-type Stirling engine
    % PROPER BETA-TYPE MODEL: Total volume depends ONLY on power piston
    % Displacer only shuttles gas between hot and cold spaces

    % In a beta-type engine:
    % - Power piston changes total gas volume
    % - Displacer moves gas between hot and cold spaces without changing total volume
    % - Total volume = f(power piston position only)

    % Total gas volume (depends ONLY on power piston position)
    % Power piston at x=0 (TDC) gives minimum volume
    % Power piston at x=2*powerCrankRadius (BDC) gives maximum volume
    totalVolume = params.totalDeadVolume + cylinderCrossSectionalArea * powerPistonPosition;

    % Regenerator volume (constant dead volume)
    regeneratorVolume = params.regeneratorVolume;

    % For Beta-type engine, we need to model how the displacer
    % shuttles gas between hot and cold spaces
    % The key is that the displacer affects the DISTRIBUTION not the TOTAL

    % Working gas volume (excluding dead volumes)
    workingGasVolume = totalVolume - params.totalDeadVolume;

    % Calculate how displacer position affects gas distribution
    % Normalize displacer position (0 at TDC, max at BDC)
    displacerMaximumStroke = 2 * displacerCrankRadius;  % Maximum stroke of displacer

    % Calculate fraction of gas in compression space based on displacer position
    % When displacerPosition = 0 (TDC), most gas is in expansion space
    % When displacerPosition = displacerMaximumStroke (BDC), most gas is in compression space
    normalizedDisplacerPosition = displacerPosition / displacerMaximumStroke;  % 0 to 1 normalized position

    % For Beta-type engine, displacer shuttles gas between hot and cold spaces
    % The displacer position affects the distribution, not the total volume
    % Use a smooth transition based on displacer position

    % Split the working volume between hot and cold spaces
    % When normalizedDisplacerPosition = 0 (TDC), more gas in expansion (hot) space
    % When normalizedDisplacerPosition = 1 (BDC), more gas in compression (cold) space

    % Use smooth sinusoidal transition for realistic gas distribution
    % This accounts for the gradual gas transfer as displacer moves
    volumeSplitFactor = 0.5 * (1 - cos(pi * normalizedDisplacerPosition));

    % Compression space (cold) volume
    compressionVolume = params.deadVolumeCold + workingGasVolume .* volumeSplitFactor;

    % Expansion space (hot) volume
    expansionVolume = params.deadVolumeHot + workingGasVolume .* (1 - volumeSplitFactor);

    % Regenerator volume stays constant
    % Note: totalVolume = compressionVolume + expansionVolume + regeneratorVolume should be satisfied

    % Alternative more physical calculation:
    % When displacer moves down, it pushes gas from cold to hot
    % compressionVolume = params.deadVolumeCold + portion of cylinder above displacer
    % expansionVolume = params.deadVolumeHot + portion of cylinder below displacer

    % Ensure conservation: totalVolume = compressionVolume + expansionVolume + regeneratorVolume
    % This is automatically satisfied by our calculation

    % Validate results
    if any(compressionVolume < 0) || any(expansionVolume < 0)
        error('Negative volume detected - check geometry parameters');
    end

    if any(totalVolume < params.totalDeadVolume)
        error('Total volume less than dead volume - check calculations');
    end

    % Calculate and display volume characteristics (only on first call)
    persistent displayed;
    if isempty(displayed)
        maximumVolume = max(totalVolume);
        minimumVolume = min(totalVolume);
        actualCompressionRatio = maximumVolume / minimumVolume;

        fprintf('Volume Analysis:\n');
        fprintf('  Max Volume: %.4f L\n', maximumVolume * 1000);
        fprintf('  Min Volume: %.4f L\n', minimumVolume * 1000);
        fprintf('  Compression Ratio: %.3f\n', actualCompressionRatio);
        fprintf('  Swept Volume (Power): %.4f L\n', params.powerSweptVolume * 1000);
        fprintf('  Swept Volume (Displacer): %.4f L\n', params.displacerSweptVolume * 1000);
        fprintf('  Dead Volume: %.4f L\n', params.totalDeadVolume * 1000);
        fprintf('\n');

        displayed = true;
    end
end