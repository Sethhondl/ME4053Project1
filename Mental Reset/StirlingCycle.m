%% Stirling Engine Cycle Problem

clear; clc; close all;

%% Functions

function pistonPosition = calculatePistonPosition(crankAngle, params, isPower)
%CALCULATEPISTONPOSITION Calculate piston position using slider-crank mechanism
%   PISTONPOSITION = CALCULATEPISTONPOSITION(CRANKANGLE, PARAMS, ISPOWER)
%   calculates the piston position for either the power piston (ISPOWER=true)
%   or the displacer (ISPOWER=false). For the displacer, the phase shift is
%   applied internally.
%
%   Inputs:
%       crankAngle - Crank angle in radians (0 = BDC, pi = TDC)
%       params     - Engine parameters structure
%       isPower    - Boolean; true for power piston, false for displacer
%
%   Output:
%       pistonPosition - Piston position in meters relative to BDC
%                       (positive values = upward from BDC)
%
%   Algorithm:
%       Uses the slider-crank kinematic equations to determine piston
%       position based on crank angle, accounting for connecting rod obliquity.
%
%   Example:
%       posP = calculatePistonPosition(pi/2, params, true);
%       posD = calculatePistonPosition(pi/2, params, false);
%
%   See also: calculateColdVolume, calculateHotVolume

    if isPower
        angle = crankAngle;
        crankLength = params.powerCrankLength;
        rodLength = params.powerRodLength;
    else
        angle = crankAngle + params.phaseShift;
        crankLength = params.displacerCrankLength;
        rodLength = params.displacerRodLength;
    end

    % Position is relative to bottom dead center (BDC)
    beta = asin(crankLength * sin(angle) / rodLength);
    pistonPosition = rodLength * cos(beta) - crankLength * cos(angle);
end

function coldVol = calculateColdVolume(crankAngle, params)
%CALCULATECOLDVOLUME Calculate cold side volume in Stirling engine
%   COLDVOL = CALCULATECOLDVOLUME(CRANKANGLE, PARAMS) calculates the cold
%   side volume of a Stirling engine at a given crank angle.
%
%   Inputs:
%       crankAngle - Crank angle in radians
%       params     - Structure containing engine parameters including:
%                   .powerCrankLength, .powerRodLength
%                   .displacerCrankLength, .displacerRodLength
%                   .phaseShift, .powerPinToPistonTop
%                   .displacerHeight, .cylinderCrossSectionalArea
%
%   Output:
%       coldVol    - Structure containing:
%                   .volume - Cold side volume in m³
%                   .height - Cold side height in m
%
%   Algorithm:
%       Calculates piston positions for both power piston and displacer,
%       then determines the cold side volume based on the space between
%       the displacer and power piston.
%
%   Example:
%       coldVol = calculateColdVolume(pi/4, engineParams);
%
%   See also: calculateHotVolume, calculatePistonPosition

    % Calculate piston positions
    powerPistonPos = calculatePistonPosition(crankAngle, params, true);
    displacerPos = calculatePistonPosition(crankAngle, params, false);
    
    % Calculate cold side height
    % Distance between displacer and power piston, minus powerPinToPistonTop, minus half displacer height
    coldVol.height = (displacerPos - powerPistonPos) - params.powerPinToPistonTop - (params.displacerHeight / 2);
    
    % Calculate cold volume
    coldVol.volume = params.cylinderCrossSectionalArea * coldVol.height;
    
    % No clamping: allow signed volume per instruction
end

function hotVol = calculateHotVolume(crankAngle, params)
%CALCULATEHOTVOLUME Calculate hot side volume in Stirling engine
%   HOTVOL = CALCULATEHOTVOLUME(CRANKANGLE, PARAMS) calculates the hot
%   side volume of a Stirling engine at a given crank angle.
%
%   Inputs:
%       crankAngle - Crank angle in radians
%       params     - Structure containing engine parameters including:
%                   .powerCrankLength, .powerRodLength
%                   .displacerCrankLength, .displacerRodLength
%                   .phaseShift, .totalCylinderHeight
%                   .displacerHeight, .cylinderCrossSectionalArea
%
%   Output:
%       hotVol     - Structure containing:
%                   .volume - Hot side volume in m³
%                   .height - Hot side height in m
%
%   Algorithm:
%       Calculates piston positions and determines the hot side volume
%       based on the space above the displacer in the cylinder.
%
%   Example:
%       hotVol = calculateHotVolume(pi/4, engineParams);
%
%   See also: calculateColdVolume, calculatePistonPosition

    % Calculate piston positions
    displacerPos = calculatePistonPosition(crankAngle, params, false);
    
    % Calculate hot side height
    % Hot height = ColdHotHeight - cold height
    hotVol.height = params.totalCylinderHeight - .5 * params.displacerHeight - displacerPos;
    
    % Calculate hot volume
    hotVol.volume = params.cylinderCrossSectionalArea * hotVol.height;
    
    % No clamping: allow signed volume per instruction
end

function schmidt = calculateSchmidtAnalysis(crankAngle, params)
%CALCULATESCHMIDTANALYSIS Calculate Schmidt analysis for Stirling engine
%   SCHMIDT = CALCULATESCHMIDTANALYSIS(CRANKANGLE, PARAMS) performs Schmidt
%   analysis to calculate instantaneous pressure and thermodynamic properties
%   of a Stirling engine at a given crank angle.
%
%   Inputs:
%       crankAngle - Crank angle in radians
%       params     - Structure containing engine parameters including:
%                   .coldTemperature, .hotTemperature, .regeneratorTemperature
%                   .regeneratorVolume, .gasConstant, .pressureAtBDC
%
%   Output:
%       schmidt    - Structure containing Schmidt analysis results:
%                   .pressure - Instantaneous pressure in Pa
%                   .totalMass - Total mass in system in kg
%                   .coldVolume - Cold side volume in m³
%                   .hotVolume - Hot side volume in m³
%                   .regeneratorVolume - Regenerator volume in m³
%
%   Algorithm:
%       Uses the Schmidt analysis method, which assumes isothermal processes
%       and ideal gas behavior. Calculates total system mass at BDC condition,
%       then uses the Schmidt equation to find instantaneous pressure.
%
%   Example:
%       schmidt = calculateSchmidtAnalysis(pi/2, engineParams);
%
%   See also: calculateColdVolume, calculateHotVolume

    % Calculate volumes at current crank angle
    coldVol = calculateColdVolume(crankAngle, params);
    hotVol = calculateHotVolume(crankAngle, params);
    
    % Extract volumes
    V_c = coldVol.volume;           % Cold side volume
    V_h = hotVol.volume;            % Hot side volume
    V_reg = params.regeneratorVolume; % Regenerator volume (constant)
    
    % Extract temperatures
    T_c = params.coldTemperature;   % Cold side temperature
    T_h = params.hotTemperature;    % Hot side temperature
    T_r = params.regeneratorTemperature; % Regenerator temperature
    
    % Extract gas constant
    R = params.gasConstant;         % Gas constant
    
    % Calculate total mass using ideal gas law at BDC condition
    % At BDC, we know the pressure and can calculate total mass
    V_comp_bdc = calculateColdVolume(0, params).volume;  % Cold volume at BDC
    V_exp_bdc = calculateHotVolume(0, params).volume;    % Hot volume at BDC
    P_bdc = params.pressureAtBDC;                        % Pressure at BDC
    
    % Apply ideal gas law at BDC condition
    denominator_bdc = V_comp_bdc/T_c + V_reg/T_r + V_exp_bdc/T_h;
    m_total = P_bdc * denominator_bdc / R;  % kg
    
    % Schmidt equation: P = m_total * R / (V_c/T_c + V_r/T_r + V_h/T_h)
    denominator = V_c/T_c + V_reg/T_r + V_h/T_h;
    P = (m_total * R) / denominator;  % Pa
    
    % Store results in structure
    schmidt.pressure = P;
    schmidt.totalMass = m_total;
    schmidt.coldVolume = V_c;
    schmidt.hotVolume = V_h;
    schmidt.regeneratorVolume = V_reg;
end

function torque = calculateTorque(crankAngle, params)
%CALCULATETORQUE Calculate engine torque at given crank angle
%   TORQUE = CALCULATETORQUE(CRANKANGLE, PARAMS) calculates the instantaneous
%   torque produced by a Stirling engine at a given crank angle.
%
%   Inputs:
%       crankAngle - Crank angle in radians
%       params     - Structure containing engine parameters including:
%                   .powerCrankLength, .powerRodLength
%                   .cylinderCrossSectionalArea, .atmosphericPressure
%
%   Output:
%       torque     - Structure containing torque results:
%                   .power - Power piston torque in N·m
%                   .displacer - Displacer torque in N·m (always 0)
%                   .total - Total engine torque in N·m
%                   .mean - Mean torque in N·m (same as total for single angle)
%
%   Algorithm:
%       Uses Schmidt analysis to get pressure, then calculates net force
%       on power piston accounting for atmospheric pressure. Applies
%       slider-crank kinematics to convert force to torque, including
%       connecting rod obliquity effects. Uses sign convention where
%       positive torque represents work-producing (expansion) phase.
%
%   Example:
%       torque = calculateTorque(pi/4, engineParams);
%
%   See also: calculateSchmidtAnalysis, calculatePistonPosition

    % Pressure from Schmidt at this crank angle (make sure its absolute)
    schmidt = calculateSchmidtAnalysis(crankAngle, params);
    P      = schmidt.pressure;

    % Geometry / constants
    r  = params.powerCrankLength;
    l  = params.powerRodLength;
    A  = params.cylinderCrossSectionalArea;   % same as cylinderArea you use elsewhere
    Patm = params.atmosphericPressure;

    % Net axial force on power piston
    Fp = (P - Patm) * A;

    % Rod obliquity
    sb = (r/l) * sin(crankAngle);             % sinβ

    cb = sqrt(1 - sb.^2);                 % cosβ > 0
    
    % Calculate torque with correct sign convention
    % Positive torque when force opposes piston motion (expansion)
    % Negative torque when force assists piston motion (compression)
    torque.power = -Fp * r * sin(crankAngle) ./ cb;

    torque.displacer = 0;                     % equal pressure both sides, zero rod area
    torque.total     = torque.power;
    torque.mean      = torque.total;          % single angle → same as total
end


function flywheel = sizeFlywheel(theta, T_total, params)
%SIZEFLYWHEEL Calculate required flywheel dimensions based on energy fluctuation
%   FLYWHEEL = SIZEFLYWHEEL(THETA, T_TOTAL, PARAMS) calculates the required
%   flywheel dimensions to achieve a specified coefficient of fluctuation
%   for a Stirling engine.
%
%   Inputs:
%       theta   - Crank angle array in radians
%       T_total - Total torque array for entire cycle in N·m
%       params  - Structure containing engine parameters including:
%                .averageRPM, .flywheelCoefficientOfFluctuation
%                .flywheelWidth, .flywheelRimThickness
%                .flywheelMaterialDensity, .maximumFlywheelDiameter
%
%   Output:
%       flywheel - Structure containing flywheel sizing results:
%                 .outerDiameter - Outer diameter in m
%                 .innerDiameter - Inner diameter in m
%                 .requiredInertia - Required moment of inertia in kg·m²
%                 .mass - Flywheel mass in kg
%                 .energyFluctuation - Energy fluctuation in J
%                 .volume - Flywheel volume in m³
%
%   Algorithm:
%       Calculates energy fluctuation from torque variation, determines
%       required moment of inertia, then uses iterative method to find
%       flywheel dimensions that provide the required inertia.
%
%   Example:
%       flywheel = sizeFlywheel(theta, totalTorque, engineParams);
%
%   See also: simulateDynamics, calculateTorque

    % Extract parameters
    omega_avg = params.averageRPM * 2*pi/60;  % Convert RPM to rad/s
    Cs = params.flywheelCoefficientOfFluctuation;
    w = params.flywheelWidth;
    t = params.flywheelRimThickness;
    rho = params.flywheelMaterialDensity;
    
    T_mean = mean(T_total);
    T_deviation = T_total - T_mean;
    
    % Calculate energy fluctuation using vectorized trapezoidal integration
    energy_variation = cumtrapz(theta, T_deviation);
    
    E_max = max(energy_variation);
    E_min = min(energy_variation);
    energy_fluctuation = E_max - E_min;
    
    % Required moment of inertia
    I_required = energy_fluctuation / (Cs * omega_avg^2);
    
    % Initial guess for outer radius using ring geometry approximation
    % This provides a better starting point than assuming a solid disk
    r_outer_guess = sqrt(I_required / (pi * rho * w * t)) + t/2;
    
    % Iterative solution to find exact flywheel dimensions
    % We need to solve: I_required = 0.5 * mass * (r_outer^2 + r_inner^2)
    % where mass = rho * volume and volume = pi * w * (r_outer^2 - r_inner^2)
    % This is non-linear, so we use iteration to converge to the solution
    
    max_iterations = params.flywheelMaxIterations;
    convergence_tolerance = params.flywheelConvergenceTolerance;  % relative error tolerance
    
    for iteration = 1:max_iterations
        % Calculate current geometry
        r_inner_current = r_outer_guess - t;
        
        % Calculate volume and mass for current radius
        flywheel_volume = pi * w * (r_outer_guess^2 - r_inner_current^2);
        flywheel_mass = rho * flywheel_volume;
        
        % Calculate actual moment of inertia for current geometry
        I_actual = 0.5 * flywheel_mass * (r_outer_guess^2 + r_inner_current^2);
        
        % Check if we have converged to the required inertia
        relative_error = abs(I_actual - I_required) / I_required;
        if relative_error < convergence_tolerance
            break;
        end
        
        % Adjust radius estimate using error ratio
        % The 1/3 power provides stable convergence for this type of problem
        error_ratio = I_required / I_actual;
        r_outer_guess = r_outer_guess * error_ratio^(1/3);
    end
    
    % Use final calculated values 
    r_outer = r_outer_guess;
    r_inner = r_outer - t;
    D_outer = 2 * r_outer;
    D_inner = 2 * r_inner;
    volume = pi * w * (r_outer^2 - r_inner^2);
    mass = rho * volume;
    
    % Store results in structure
    flywheel.outerDiameter = D_outer;
    flywheel.innerDiameter = D_inner;
    flywheel.requiredInertia = I_required;
    flywheel.mass = mass;
    flywheel.energyFluctuation = energy_fluctuation;
    flywheel.volume = volume;
    
end

function dynamics = simulateDynamics(theta, T_total, I_flywheel, params)
%SIMULATEDYNAMICS Simulate angular velocity variation with flywheel
%   DYNAMICS = SIMULATEDYNAMICS(THETA, T_TOTAL, I_FLYWHEEL, PARAMS) simulates
%   the angular velocity variation of a Stirling engine with a flywheel
%   over one complete cycle.
%
%   Inputs:
%       theta      - Crank angle array in radians
%       T_total    - Total torque array for entire cycle in N·m
%       I_flywheel - Flywheel moment of inertia in kg·m²
%       params     - Structure containing engine parameters including:
%                   .averageRPM, .flywheelCoefficientOfFluctuation
%
%   Output:
%       dynamics   - Structure containing dynamics simulation results:
%                   .angularVelocity - Angular velocity in rad/s
%                   .angularAcceleration - Angular acceleration in rad/s²
%                   .rpm - Angular velocity in RPM
%                   .coefficientOfFluctuation - Actual coefficient of fluctuation
%                   .netTorque - Net torque (engine - load) in N·m
%                   .loadTorque - Load torque in N·m
%
%   Algorithm:
%       Uses work-energy theorem to calculate velocity variation based on
%       net torque (engine torque minus load torque). Load torque equals
%       mean engine torque for steady-state operation.
%
%   Example:
%       dynamics = simulateDynamics(theta, totalTorque, flywheelInertia, engineParams);
%
%   See also: sizeFlywheel, calculateTorque

    % Extract parameters
    omega_target = params.averageRPM * 2*pi/60;  % Convert RPM to rad/s
    
    % Calculate load and net torques
    % Load torque equals mean engine torque (steady-state condition)
    T_load = mean(T_total);
    T_net = T_total - T_load;  % Net torque available for acceleration/deceleration
    
    % Calculate angular acceleration from net torque
    angular_acceleration = T_net / I_flywheel;
    
    % Energy-based velocity calculation
    % Use work-energy theorem: dW = 0.5 * I * (omega_f^2 - omega_i^2)
    % where dW is work done by net torque over angle increment
    
    % Calculate cumulative work done by net torque using vectorized integration
    cumulative_work = cumtrapz(theta, T_net);
    
    % Apply work-energy theorem to find velocity at each point
    % Starting from omega_target, work changes the kinetic energy
    velocity_squared = omega_target^2 + 2 * cumulative_work / I_flywheel;
    
    % No minimum-speed floor: compute directly from energy balance
    angular_velocity = sqrt(velocity_squared);
    
    % Normalize to maintain correct average speed
    % This corrects for any drift in the energy-based calculation
    omega_actual_avg = mean(angular_velocity);
    angular_velocity = angular_velocity * (omega_target / omega_actual_avg);
    
    % Convert to RPM for convenience
    rpm = angular_velocity * 60 / (2*pi);
    
    % Calculate actual coefficient of fluctuation
    omega_max = max(angular_velocity);
    omega_min = min(angular_velocity);
    omega_mean = mean(angular_velocity);
    coefficient_of_fluctuation = (omega_max - omega_min) / omega_mean;
    
    % Store results in structure
    dynamics.angularVelocity = angular_velocity;
    dynamics.angularAcceleration = angular_acceleration;
    dynamics.rpm = rpm;
    dynamics.coefficientOfFluctuation = coefficient_of_fluctuation;
    dynamics.netTorque = T_net;
    dynamics.loadTorque = T_load;
end



function optimization = optimizePhaseShift(theta, params)
%OPTIMIZEPHASESHIFT Find phase shift that maximizes power output
%   OPTIMIZATION = OPTIMIZEPHASESHIFT(THETA, PARAMS) performs a three-stage
%   search over phase shift to maximize power. It scans broadly, then narrows
%   around the best candidate with increasingly fine resolution down to 0.01°.
%
%   Inputs:
%       theta  - Crank angle array in radians
%       params - Engine parameters structure (uses .averageRPM)
%
%   Output (structure):
%       optimization.phaseGridCoarse      - Tested coarse grid (rad)
%       optimization.meanTorqueCoarse     - Mean torque per coarse phase (N·m)
%       optimization.powerCoarse          - Power per coarse phase (W)
%       optimization.phaseGridMedium      - Tested medium grid (rad)
%       optimization.meanTorqueMedium     - Mean torque per medium phase (N·m)
%       optimization.powerMedium          - Power per medium phase (W)
%       optimization.phaseGridFine        - Tested fine grid (rad)
%       optimization.meanTorqueFine       - Mean torque per fine phase (N·m)
%       optimization.powerFine            - Power per fine phase (W)
%       optimization.bestPhaseShift       - Best phase (rad, from fine scan)
%       optimization.bestPower            - Max power (W, from fine scan)
%       optimization.bestMeanTorque       - Mean torque at best phase (N·m)
%
%   Notes:
%       Power is computed as meanTorque * omega_avg, which is equivalent to
%       integrating torque over angle and multiplying by rotational speed.

    omega_avg = params.averageRPM * 2*pi/60;  % rad/s

    % ---------- Stage 1: Coarse scan (broad range) ----------
    phaseGridCoarse = deg2rad(30):deg2rad(2):deg2rad(150); % 30° to 150° in 2° steps
    [meanTorqueCoarse, powerCoarse] = evaluateGrid(theta, params, phaseGridCoarse, omega_avg);

    [~, idxCoarse] = max(powerCoarse);
    center1 = phaseGridCoarse(idxCoarse);

    % ---------- Stage 2: Medium scan (narrow window) ----------
    window2 = deg2rad(6);                % ±6° around coarse best
    step2   = deg2rad(0.1);              % 0.1° resolution
    phaseGridMedium = (center1 - window2):step2:(center1 + window2);
    [meanTorqueMedium, powerMedium] = evaluateGrid(theta, params, phaseGridMedium, omega_avg);

    [~, idxMedium] = max(powerMedium);
    center2 = phaseGridMedium(idxMedium);

    % ---------- Stage 3: Fine scan (very narrow, 0.01°) ----------
    window3 = deg2rad(0.5);              % ±0.5° around medium best
    step3   = deg2rad(0.01);             % 0.01° resolution
    phaseGridFine = (center2 - window3):step3:(center2 + window3);
    [meanTorqueFine, powerFine] = evaluateGrid(theta, params, phaseGridFine, omega_avg);

    [bestPower, idxFine] = max(powerFine);
    bestPhaseShift = phaseGridFine(idxFine);
    bestMeanTorque = meanTorqueFine(idxFine);

    % Package results
    optimization.phaseGridCoarse    = phaseGridCoarse;
    optimization.meanTorqueCoarse   = meanTorqueCoarse;
    optimization.powerCoarse        = powerCoarse;
    optimization.phaseGridMedium    = phaseGridMedium;
    optimization.meanTorqueMedium   = meanTorqueMedium;
    optimization.powerMedium        = powerMedium;
    optimization.phaseGridFine      = phaseGridFine;
    optimization.meanTorqueFine     = meanTorqueFine;
    optimization.powerFine          = powerFine;
    optimization.bestPhaseShift     = bestPhaseShift;
    optimization.bestPower          = bestPower;
    optimization.bestMeanTorque     = bestMeanTorque;

    function [meanTorqueArr, powerArr] = evaluateGrid(thetaLoc, paramsLoc, gridRad, omegaAvg)
        meanTorqueArr = zeros(size(gridRad));
        powerArr = zeros(size(gridRad));
        for kk = 1:numel(gridRad)
            paramsLoc.phaseShift = gridRad(kk);
            T_total_loc = zeros(size(thetaLoc));
            for ii = 1:length(thetaLoc)
                tqLoc = calculateTorque(thetaLoc(ii), paramsLoc);
                T_total_loc(ii) = tqLoc.total;
            end
            mT = mean(T_total_loc);
            meanTorqueArr(kk) = mT;
            powerArr(kk) = mT * omegaAvg;
        end
    end
end



%% Prescribed Parameters

% ============== GEOMETRY PARAMETERS ==============

% Power Piston
params.powerCrankLength = 0.025;              % m - Crank radius
params.powerRodLength = 0.075;                % m - Connecting rod length
params.powerPinToPistonTop = 0.005;           % m - Pin to piston distance

% Displacer
params.displacerCrankLength = 0.020;          % m - Crank radius
params.displacerRodLength = 0.140;            % m - Connecting rod length
params.displacerVolume = 4e-5;                % m³ - Displacer volume

% Cylinder
params.cylinderBore = 0.050;                  % m - Bore diameter



% ============== OPERATING CONDITIONS ==============

% Kinematics
params.phaseShift = pi/2;                     % rad - Phase angle (90°)
params.averageRPM = 650;                      % RPM - Operating speed

% Thermodynamics
params.hotTemperature = 900;                  % K - Hot space temperature
params.coldTemperature = 300;                 % K - Cold space temperature


% Pressure
params.pressureAtBDC = 500e3;                 % Pa - Pressure at bottom dead center
params.atmosphericPressure = 101.3e3;         % Pa - Atmospheric pressure


% ============== DEAD VOLUMES ==============

params.compressionRatio = 1.7;                % Given compression ratio
params.regeneratorVolume = 2e-5;              % m³ - Given regenerator volume


% ============== WORKING FLUID (AIR) ==============

params.gasConstant = 287;                     % J/(kg·K) - Specific gas constant
params.gasGamma = 1.4;                        % - Heat capacity ratio
params.gasName = 'Air';


% ============== FLYWHEEL SPECIFICATIONS ==============

params.flywheelWidth = 0.025;                 % m - Flywheel width
params.flywheelRimThickness = 0.050;          % m - Rim thickness
params.flywheelMaterialDensity = 7750;        % kg/m³ - 304 Stainless Steel
params.flywheelCoefficientOfFluctuation = 0.003;  % - Target Cs value

% ============== SIMULATION PARAMETERS ==============

params.simulationPointsPerCycle = 360;        % Points per cycle
params.simulationCycles = 3;                  % Number of cycles
params.simulationTolerance = 1e-6;            % Convergence tolerance
params.maximumFlywheelDiameter = 2.0;         % m - Maximum allowable flywheel diameter
params.flywheelMaxIterations = 20;            % Iterations for flywheel sizing
params.flywheelConvergenceTolerance = 1e-3;   % Relative error tolerance for inertia match


%% Assisting Calculations

% Cylinder cross-sectional area
params.cylinderCrossSectionalArea = pi/4*(params.cylinderBore)^2;

% Displacer GEOMETRY
params.displacerHeight = params.displacerVolume / params.cylinderCrossSectionalArea;

% Calculate power piston positions at BDC and TDC
params.powerPistonPosBDC = calculatePistonPosition(0, params, true);  % 0 radians = BDC
params.powerPistonPosTDC = calculatePistonPosition(pi, params, true);  % π radians = TDC

% Calculate swept volume (volume displaced by power piston)
params.powerSweptVolume = params.cylinderCrossSectionalArea * (params.powerPistonPosTDC - params.powerPistonPosBDC);

% Calculate total volume at BDC (9.23 Pg 2) (Also V_MAX)
params.totalVolumeBDC = params.regeneratorVolume - params.displacerVolume + ( params.compressionRatio * params.powerSweptVolume ) / ( params.compressionRatio - 1 ) ;

%Total cylinder height (crank pin to cylinder roof)
params.ColdHotHeight = params.totalVolumeBDC / params.cylinderCrossSectionalArea;
params.totalCylinderHeight = params.ColdHotHeight + params.displacerHeight + params.powerPinToPistonTop + params.powerRodLength -params.powerCrankLength;

% Regenerator temperature
params.regeneratorTemperature = (params.hotTemperature + params.coldTemperature) / 2;  % K












%% Test Animation - Stirling Engine Cycle Analysis

% Create crank angle array for one complete cycle
theta = linspace(0, 2*pi, params.simulationPointsPerCycle);

% Initialize cycle data structure
cycleData.powerPistonPos = zeros(size(theta));
cycleData.displacerPos = zeros(size(theta));
cycleData.totalVolume = zeros(size(theta));
cycleData.hotVolume = zeros(size(theta));
cycleData.coldVolume = zeros(size(theta));
cycleData.regeneratorVolume = params.regeneratorVolume * ones(size(theta));
cycleData.pressure = zeros(size(theta));
cycleData.totalTorque = zeros(size(theta));
cycleData.powerTorque = zeros(size(theta));

% Calculate all data for each crank angle
for i = 1:length(theta)
    % Calculate piston positions
    cycleData.powerPistonPos(i) = calculatePistonPosition(theta(i), params, true);
    cycleData.displacerPos(i) = calculatePistonPosition(theta(i), params, false);
    
    % Calculate volumes
    coldVol = calculateColdVolume(theta(i), params);
    hotVol = calculateHotVolume(theta(i), params);
    
    cycleData.coldVolume(i) = coldVol.volume;
    cycleData.hotVolume(i) = hotVol.volume;
    cycleData.totalVolume(i) = cycleData.coldVolume(i) + cycleData.hotVolume(i) + cycleData.regeneratorVolume(i);
    
    % Calculate Schmidt analysis
    schmidt = calculateSchmidtAnalysis(theta(i), params);
    cycleData.pressure(i) = schmidt.pressure;
    
    % Calculate torque
    torque = calculateTorque(theta(i), params);
    cycleData.totalTorque(i) = torque.total;
    cycleData.powerTorque(i) = torque.power;
end

% Calculate flywheel sizing and dynamics
flywheel = sizeFlywheel(theta, cycleData.totalTorque, params);
dynamics = simulateDynamics(theta, cycleData.totalTorque, flywheel.requiredInertia, params);

% Calculate cycle statistics and store in results structure
results.meanPressure = mean(cycleData.pressure);
results.maxPressure = max(cycleData.pressure);
results.minPressure = min(cycleData.pressure);
results.meanTorque = mean(cycleData.totalTorque);
results.meanAngularVelocity = mean(dynamics.rpm);

% Phase optimization (scan phase shift to maximize power)
optimization = optimizePhaseShift(theta, params);

% Create figures for analysis
figure('Name', 'Piston Positions vs Crank Angle', 'Position', [50, 50, 900, 600]);
plot(theta*180/pi, cycleData.powerPistonPos*1000, 'b-', 'LineWidth', 2, 'DisplayName', 'Power Piston');
hold on;
plot(theta*180/pi, cycleData.displacerPos*1000, 'r-', 'LineWidth', 2, 'DisplayName', 'Displacer');
xlabel('Crank Angle (degrees)');
ylabel('Position (mm)');
title('Piston Positions vs Crank Angle');
legend('Location', 'best');
grid on;

% Phase Optimization Results - separate figures (use fine grid for display)
phaseDeg = optimization.phaseGridFine * 180/pi;
bestPhaseDeg = optimization.bestPhaseShift * 180/pi;

figure('Name', 'Power vs Phase Shift', 'Position', [230, 230, 900, 600]);
plot(phaseDeg, optimization.powerFine, 'k-', 'LineWidth', 2, 'DisplayName', 'Power');
hold on;
plot(bestPhaseDeg, optimization.bestPower, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Best');
xlabel('Phase Shift (degrees)');
ylabel('Power (W)');
title('Power vs Phase Shift');
legend('Location', 'best');
grid on;

figure('Name', 'Mean Torque vs Phase Shift', 'Position', [260, 260, 900, 600]);
plot(phaseDeg, optimization.meanTorqueFine, 'b-', 'LineWidth', 2, 'DisplayName', 'Mean Torque');
hold on;
plot(bestPhaseDeg, optimization.bestMeanTorque, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Best');
xlabel('Phase Shift (degrees)');
ylabel('Mean Torque (N·m)');
title('Mean Torque vs Phase Shift');
legend('Location', 'best');
grid on;

figure('Name', 'Volumes vs Crank Angle', 'Position', [80, 80, 900, 600]);
plot(theta*180/pi, cycleData.totalVolume*1e6, 'k-', 'LineWidth', 2, 'DisplayName', 'Total Volume');
hold on;
plot(theta*180/pi, cycleData.hotVolume*1e6, 'r-', 'LineWidth', 2, 'DisplayName', 'Hot Volume');
plot(theta*180/pi, cycleData.coldVolume*1e6, 'b-', 'LineWidth', 2, 'DisplayName', 'Cold Volume');
plot(theta*180/pi, cycleData.regeneratorVolume*1e6, 'g-', 'LineWidth', 2, 'DisplayName', 'Regenerator Volume');
xlabel('Crank Angle (degrees)');
ylabel('Volume (cm³)');
title('Volumes vs Crank Angle');
legend('Location', 'best');
grid on;

figure('Name', 'Pressure vs Crank Angle', 'Position', [110, 110, 900, 600]);
plot(theta*180/pi, cycleData.pressure/1000, 'k-', 'LineWidth', 2, 'DisplayName', 'Instantaneous Pressure');
hold on;
plot(theta*180/pi, results.meanPressure/1000*ones(size(theta)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Mean Pressure');
plot(theta*180/pi, results.maxPressure/1000*ones(size(theta)), 'g--', 'LineWidth', 1, 'DisplayName', 'Max Pressure');
plot(theta*180/pi, results.minPressure/1000*ones(size(theta)), 'b--', 'LineWidth', 1, 'DisplayName', 'Min Pressure');
xlabel('Crank Angle (degrees)');
ylabel('Pressure (kPa)');
title('Pressure vs Crank Angle');
legend('Location', 'best');
grid on;

% Compute specific volume (total volume per total mass)
m_total_system = calculateSchmidtAnalysis(0, params).totalMass;  % kg (constant over cycle)
specificVolume = cycleData.totalVolume / m_total_system;          % m^3/kg

figure('Name', 'P-v Diagram (Schmidt Analysis)', 'Position', [140, 140, 900, 600]);
plot(specificVolume*1e6, cycleData.pressure/1000, 'k-', 'LineWidth', 2, 'DisplayName', 'Cycle');
hold on;
% Add constant pressure lines
plot([min(specificVolume)*1e6, max(specificVolume)*1e6], [results.meanPressure/1000, results.meanPressure/1000], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Mean Pressure');
plot([min(specificVolume)*1e6, max(specificVolume)*1e6], [results.maxPressure/1000, results.maxPressure/1000], 'g--', 'LineWidth', 1, 'DisplayName', 'Max Pressure');
plot([min(specificVolume)*1e6, max(specificVolume)*1e6], [results.minPressure/1000, results.minPressure/1000], 'b--', 'LineWidth', 1, 'DisplayName', 'Min Pressure');
xlabel('Specific Volume (cm^3/kg)');
ylabel('Pressure (kPa)');
title('P-v Diagram (Schmidt Analysis)');
legend('Location', 'best');
grid on;
axis tight;

figure('Name', 'Torque vs Crank Angle', 'Position', [170, 170, 900, 600]);
plot(theta*180/pi, cycleData.totalTorque, 'k-', 'LineWidth', 2, 'DisplayName', 'Total Torque');
hold on;
plot(theta*180/pi, cycleData.powerTorque, 'b-', 'LineWidth', 2, 'DisplayName', 'Power Piston Torque');
plot(theta*180/pi, results.meanTorque*ones(size(theta)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Mean Torque');
xlabel('Crank Angle (degrees)');
ylabel('Torque (N·m)');
title('Torque vs Crank Angle');
legend('Location', 'best');
grid on;

figure('Name', 'Angular Velocity vs Crank Angle', 'Position', [200, 200, 900, 600]);
plot(theta*180/pi, dynamics.rpm, 'k-', 'LineWidth', 2, 'DisplayName', 'Angular Velocity');
hold on;
plot(theta*180/pi, results.meanAngularVelocity*ones(size(theta)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Mean Velocity');
plot(theta*180/pi, params.averageRPM*ones(size(theta)), 'g--', 'LineWidth', 1, 'DisplayName', 'Target Velocity');
xlabel('Crank Angle (degrees)');
ylabel('Angular Velocity (RPM)');
title('Angular Velocity vs Crank Angle');
legend('Location', 'best');
grid on;

% Display comprehensive analysis results
fprintf('\n=== STIRLING ENGINE CYCLE ANALYSIS RESULTS ===\n');
fprintf('Engine Configuration:\n');
fprintf('  Power Piston Crank Length: %.1f mm\n', params.powerCrankLength*1000);
fprintf('  Power Piston Rod Length: %.1f mm\n', params.powerRodLength*1000);
fprintf('  Displacer Crank Length: %.1f mm\n', params.displacerCrankLength*1000);
fprintf('  Displacer Rod Length: %.1f mm\n', params.displacerRodLength*1000);
fprintf('  Phase Shift: %.1f degrees\n', params.phaseShift*180/pi);
fprintf('  Cylinder Bore: %.1f mm\n', params.cylinderBore*1000);

fprintf('\nThermodynamic Conditions:\n');
fprintf('  Hot Temperature: %.0f K (%.0f°C)\n', params.hotTemperature, params.hotTemperature-273.15);
fprintf('  Cold Temperature: %.0f K (%.0f°C)\n', params.coldTemperature, params.coldTemperature-273.15);
fprintf('  Regenerator Temperature: %.0f K (%.0f°C)\n', params.regeneratorTemperature, params.regeneratorTemperature-273.15);
fprintf('  Pressure at BDC: %.1f kPa\n', params.pressureAtBDC/1000);

fprintf('\nPiston Motion Analysis:\n');
fprintf('  Power Piston Stroke: %.2f mm\n', (max(cycleData.powerPistonPos) - min(cycleData.powerPistonPos))*1000);
fprintf('  Displacer Stroke: %.2f mm\n', (max(cycleData.displacerPos) - min(cycleData.displacerPos))*1000);
fprintf('  Power Piston Max Position: %.2f mm\n', max(cycleData.powerPistonPos)*1000);
fprintf('  Displacer Max Position: %.2f mm\n', max(cycleData.displacerPos)*1000);

fprintf('\nVolume Analysis:\n');
fprintf('  Total Volume Range: %.2f - %.2f cm³\n', min(cycleData.totalVolume)*1e6, max(cycleData.totalVolume)*1e6);
fprintf('  Compression Ratio: %.2f\n', max(cycleData.totalVolume)/min(cycleData.totalVolume));
fprintf('  Hot Volume Range: %.2f - %.2f cm³\n', min(cycleData.hotVolume)*1e6, max(cycleData.hotVolume)*1e6);
fprintf('  Cold Volume Range: %.2f - %.2f cm³\n', min(cycleData.coldVolume)*1e6, max(cycleData.coldVolume)*1e6);
fprintf('  Regenerator Volume: %.2f cm³ (constant)\n', cycleData.regeneratorVolume(1)*1e6);

fprintf('\nPressure Analysis (Schmidt):\n');
fprintf('  Pressure Range: %.2f - %.2f kPa\n', min(cycleData.pressure)/1000, max(cycleData.pressure)/1000);
fprintf('  Mean Pressure: %.2f kPa\n', results.meanPressure/1000);
fprintf('  Pressure Ratio: %.2f\n', max(cycleData.pressure)/min(cycleData.pressure));

fprintf('\nTorque Analysis:\n');
fprintf('  Total Torque Range: %.3f - %.3f N·m\n', min(cycleData.totalTorque), max(cycleData.totalTorque));
fprintf('  Mean Total Torque: %.3f N·m\n', results.meanTorque);
fprintf('  Power Piston Torque Range: %.3f - %.3f N·m\n', min(cycleData.powerTorque), max(cycleData.powerTorque));
fprintf('  Mean Power Piston Torque: %.3f N·m\n', mean(cycleData.powerTorque));

fprintf('\nFlywheel Analysis:\n');
fprintf('  Required Moment of Inertia: %.4f kg·m²\n', flywheel.requiredInertia);
fprintf('  Flywheel Outer Diameter: %.3f m (%.1f mm)\n', flywheel.outerDiameter, flywheel.outerDiameter*1000);
fprintf('  Flywheel Inner Diameter: %.3f m (%.1f mm)\n', flywheel.innerDiameter, flywheel.innerDiameter*1000);
fprintf('  Flywheel Mass: %.2f kg\n', flywheel.mass);
fprintf('  Flywheel Volume: %.6f m³ (%.1f cm³)\n', flywheel.volume, flywheel.volume*1e6);
fprintf('  Energy Fluctuation: %.3f J\n', flywheel.energyFluctuation);

fprintf('\nDynamics Analysis:\n');
fprintf('  Angular Velocity Range: %.1f - %.1f RPM\n', min(dynamics.rpm), max(dynamics.rpm));
fprintf('  Mean Angular Velocity: %.1f RPM\n', results.meanAngularVelocity);
fprintf('  Target Angular Velocity: %.1f RPM\n', params.averageRPM);
fprintf('  Actual Coefficient of Fluctuation: %.4f\n', dynamics.coefficientOfFluctuation);
fprintf('  Target Coefficient of Fluctuation: %.4f\n', params.flywheelCoefficientOfFluctuation);
fprintf('  Net Torque Range: %.3f - %.3f N·m\n', min(dynamics.netTorque), max(dynamics.netTorque));
fprintf('  Load Torque: %.3f N·m\n', dynamics.loadTorque);

fprintf('\nPerformance Summary:\n');
fprintf('  Engine Power Output: %.2f W (estimated from mean torque)\n', results.meanTorque * params.averageRPM * 2*pi/60);
fprintf('  Flywheel Effectiveness: %.1f%% (target Cs: %.4f, actual Cs: %.4f)\n', ...
        (1 - dynamics.coefficientOfFluctuation/params.flywheelCoefficientOfFluctuation)*100, ...
        params.flywheelCoefficientOfFluctuation, dynamics.coefficientOfFluctuation);
fprintf('  Cycle Completeness: %.1f%% (pressure returns to within 1%% of start)\n', ...
        (1 - abs(cycleData.pressure(end) - cycleData.pressure(1))/cycleData.pressure(1))*100);

fprintf('\nPhase Optimization:\n');
fprintf('  Best Phase Shift: %.1f degrees\n', optimization.bestPhaseShift*180/pi);
fprintf('  Max Power: %.2f W\n', optimization.bestPower);
fprintf('  Mean Torque at Best Phase: %.3f N·m\n', optimization.bestMeanTorque);

fprintf('\n===============================================\n');
fprintf('\n===============================================\n');