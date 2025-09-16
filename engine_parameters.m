function params = engine_parameters()
    % STIRLING ENGINE PARAMETERS
    % Configurable parameters for Beta-type Stirling Engine Analysis
    % Large-scale power generation system (1-10 kW range)
    % Working Fluid: Air
    % Flywheel Material: Steel
    
    %% Power Piston Parameters
    params.power.crank_length = 0.060;        % m (60 mm) - optimized for 1-10kW
    params.power.rod_length = 0.240;          % m (240 mm)

    %% Displacer Parameters
    params.disp.crank_length = 0.075;         % m (75 mm) - 1.25x power stroke
    params.disp.rod_length = 0.300;           % m (300 mm)
    params.disp.volume = 0.002;               % m^3 (2000 cm^3)

    %% Cylinder Parameters
    params.cylinder.bore = 0.180;             % m (180 mm diameter) - sized for 1-10kW power
    params.cylinder.area = pi * (params.cylinder.bore/2)^2;  % m^2
    
    %% Operating Parameters
    params.phase_shift = 90 * pi/180;         % radians (90 degrees)
    params.compression_ratio = 3.5;           % dimensionless
    
    %% Temperature Conditions
    params.T_hot = 600;                       % K (hot temperature) - reduced for smaller engine
    params.T_cold = 350;                      % K (cold temperature)
    params.T_regenerator = (params.T_hot + params.T_cold) / 2;  % K (average)
    
    %% Pressure Conditions
    params.P_bdc = 1.5e6;                     % Pa (1.5 MPa at bottom dead center) - for target power
    params.P_atm = 101325;                    % Pa (atmospheric pressure)
    
    %% Dead Volumes
    params.V_dead_hot = 0.00005;              % m^3 (50 cm^3 hot space dead volume)
    params.V_dead_cold = 0.00005;             % m^3 (50 cm^3 cold space dead volume)
    params.V_regenerator = 0.0001;            % m^3 (100 cm^3 regenerator dead volume)
    params.V_dead_total = params.V_dead_hot + params.V_dead_cold + params.V_regenerator;
    
    %% Working Fluid Properties (Air)
    params.gas.R = 287;                       % J/(kg*K) - specific gas constant for air
    params.gas.gamma = 1.4;                   % heat capacity ratio
    params.gas.cv = params.gas.R / (params.gas.gamma - 1);  % J/(kg*K)
    params.gas.cp = params.gas.gamma * params.gas.cv;       % J/(kg*K)
    params.gas.name = 'Air';
    
    %% Flywheel Parameters
    params.flywheel.width = 0.100;            % m (100 mm)
    params.flywheel.thickness = 0.040;        % m (40 mm rim thickness)
    params.flywheel.material_density = 7850;  % kg/m^3 (steel)
    params.flywheel.coeff_fluctuation = 0.04; % dimensionless (4%)
    
    %% Operating Speed
    params.rpm_avg = 500;                     % RPM (average rotational speed)
    params.omega_avg = params.rpm_avg * 2*pi/60;  % rad/s
    params.frequency = params.rpm_avg / 60;   % Hz
    
    %% Simulation Parameters
    params.sim.n_points = 360;                % number of points per cycle
    params.sim.n_cycles = 3;                  % number of cycles to simulate
    params.sim.tolerance = 1e-6;              % convergence tolerance
    
    %% Validation Limits (for error checking)
    params.limits.min_efficiency = 0.0;       % minimum allowable efficiency
    params.limits.max_efficiency = 1 - params.T_cold/params.T_hot;  % Carnot limit
    params.limits.min_power = 0;              % W
    params.limits.max_power = 50000;          % W (50 kW max)
    params.limits.max_flywheel_diameter = 2.0; % m
    params.limits.min_pressure = 0;           % Pa (must be positive)
    
    %% Calculated Derived Parameters
    % Swept volumes
    params.V_swept_power = params.cylinder.area * 2 * params.power.crank_length;  % m^3
    params.V_swept_disp = params.cylinder.area * 2 * params.disp.crank_length;    % m^3
    
    % For beta-type engine, compression ratio is based on power piston swept volume
    % Maximum volume when power piston at BDC
    params.V_max = params.V_swept_power + params.V_dead_total;
    % Minimum volume when power piston at TDC
    params.V_min = params.V_dead_total;

    % Check compression ratio
    calculated_CR = params.V_max / params.V_min;
    if abs(calculated_CR - params.compression_ratio) > 0.1
        % Use the calculated value
        params.compression_ratio = calculated_CR;
    end
    
    %% Display Configuration Summary
    fprintf('\n=== STIRLING ENGINE CONFIGURATION ===\n');
    fprintf('Engine Type: Beta-type Stirling\n');
    fprintf('Working Fluid: %s\n', params.gas.name);
    fprintf('Cylinder Bore: %.0f mm\n', params.cylinder.bore * 1000);
    fprintf('Power Stroke: %.0f mm\n', params.power.crank_length * 2 * 1000);
    fprintf('Displacer Stroke: %.0f mm\n', params.disp.crank_length * 2 * 1000);
    fprintf('Phase Shift: %.0f degrees\n', params.phase_shift * 180/pi);
    fprintf('Compression Ratio: %.2f\n', params.compression_ratio);
    fprintf('Operating Speed: %.0f RPM\n', params.rpm_avg);
    fprintf('Temperature Range: %.0f K to %.0f K\n', params.T_cold, params.T_hot);
    fprintf('Carnot Efficiency Limit: %.1f%%\n', params.limits.max_efficiency * 100);
    fprintf('=====================================\n\n');
end