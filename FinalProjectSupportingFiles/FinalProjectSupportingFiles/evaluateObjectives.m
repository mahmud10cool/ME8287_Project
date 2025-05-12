function evaluatedChrom = evaluateObjectives(parameter)
    % EVALUATEDESIGN evaluates a design (candidate design) and computes
    % the objectives. Any constraints/constants if needed, can be defined local to this
    % function. FEA runs if any, will be performed here.

    %% Initialization
    %clear;
    %close all;
    addpath('C:\femm42\mfiles'); % This is the FEMM path. You may need to modify this to the installation directory of FEMM on your PC.
    %clc;
        
    %% Stator Dimensions
    dimensions.stator.dso = 2;                   %[mm]
    dimensions.stator.dsp = 4;                   %[mm]
    dimensions.stator.dst = parameter(4);        %[mm]
    dimensions.stator.dsy = parameter(3);        %[mm]
    dimensions.stator.alpha_st = parameter(6);   %[degree]
    dimensions.stator.wst = parameter(5);        %[mm]
    dimensions.outerRadius = 180;                %[mm]
    dimensions.delta = parameter(2);             %[mm] Airgap length
    
    %% Define the number of pole pairs
    p = 2;

    %% Rotor Dimensions 
    dimensions.rotor.dm=parameter(1);                    %[mm] 
    rsi = dimensions.outerRadius - dimensions.stator.dsy ...
        - dimensions.stator.dst - dimensions.stator.dsp; %[mm]
    dimensions.rotor.dri = rsi - dimensions.delta ...
        - dimensions.rotor.dm;                           %[mm]
    dimensions.rotor.dmp=dimensions.rotor.dm;            %[mm]
    dimensions.rotor.alpha_m = 180/p;                    %[deg]
        
    %% Materials 
    materials.coil.name = '18 AWG';
    materials.coil.conductivity = 5.77e7;   % Copper conductivity in S/m
    materials.coil.AWG = 18;
    materials.iron.name = 'M19_29Ga';
    materials.iron.ch = 0.0186;             % Hysteresis loss coefficient for M19-29Ga in (Watt/(kg * T^2 * Hz^2)
    materials.iron.ce = 6.887e-5;           % Eddy loss Current coefficientfor M19-29Ga in (Watts/(kg * T^2 * Hz)
    materials.iron.sf = 0.91;               % Lamination stacking factor (nondimensional)
    materials.iron.density = 7650;          % Mass density of M19 29Ga in [kg/m^3]
    materials.magnet.name = 'N40';
    materials.magnet.conductivity = 5.55e5; % Magnet conductivity in S/m
    
    steelCost = 14.03e3;                    % Cost of steel [$/m^3]
    cuCost = 0.06;                          % Cost of AWG 18 [$/m^3]
    magnetCost = 708.5e3;                   % Cost of magnet [$/m^3]
    k_ov = 1.8; 
        
    %% Define the winding structure and add parameters
    Q = 12;

    % Slot angle [rad]
    alpha_c = 2*pi/Q; 
    Kcu = 0.4;

    % Wire diameter [m]
    dw=0.324861*0.0254*exp(-0.115942*materials.coil.AWG); 
    
    % Cross sectional area [mm^2]
    Sc = 0.25*pi*dw^2*1e6; 
    Su = alpha_c*(dimensions.outerRadius-dimensions.stator.dsy)*dimensions.stator.dst - dimensions.stator.wst*dimensions.stator.dst - 0.5*alpha_c*dimensions.stator.dst^2;
    
    % Number of turns
    zQ = round(Kcu/Sc/2*Su);                    
    
    % Winding
    winding.layers = 2; 
    winding.topSlots.circuit = [{'U'}, {'W'}, {'V'}, {'U'}, {'W'}, {'V'}, {'U'}, {'W'}, {'V'}, {'U'}, {'W'}, {'V'}];
    winding.topSlots.zQ = zQ.*[1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1]; % Positive number for conductors that go into the page
    winding.bottomSlots.zQ = zQ.*[1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1];
    winding.bottomSlots.circuit = [{'U'}, {'W'}, {'V'}, {'U'}, {'W'}, {'V'}, {'U'}, {'W'}, {'V'}, {'U'}, {'W'}, {'V'}];
    winding.y = 3;
    
    %% Add parameters to the settings structure
    settings.RPM = 10000; % Rated Speed in RPM
    settings.phi = 0;    % Angle in [deg]
    settings.iHat = 5*Sc*sqrt(2);  % Peak current [A]
    settings.lowestHarmonic = 1; %Lowest Harmonic (n) for which winding factor |Kw,n| is non-zero; The max value it can take is p
    settings.steps = 14; % Number of steps the design must be evaluated at; Keep this an even number;
    
    %% Rated torque
    P = 50e3; % [W]
    ratedTorque = P/(settings.RPM*(2*pi)/60);
    
    %% Analysis
    [C, Ceq] = evaluateConstraints(parameter); % check if constaints are violated
    
    if C(1) > 0 || C(2) > 0 || C(3) > 0 || C(4) > 0 || Ceq ~= 0
        length = inf;
        designEval.torque.ripple = inf;
        designEval.efficiency = -inf;
    else
        [designEval, length] = evaluateDesign(materials, dimensions, p, winding, settings, ratedTorque);
        fprintf('The avg. torque is %1.3f [Nm]; the expected value is 47.7465[Nm]\n', designEval.torque.average);
        fprintf('The Efficiency is %1.3f percent; the expected value is 90.00 percent\n', designEval.efficiency);
        fprintf('The torque ripple is %1.3f; the expected value is 0.50\n', designEval.torque.ripple);
    end
    
       
    %% Compute Objectives.
   
    % Axial length from evaluateDesign
    l = length; 

    % Rotor steel volume [m^3]
    Vrs = pi*(rsi - dimensions.delta - dimensions.rotor.dm)^2*l/1e9; 
    
    % Stator steel volume [m^3]
    Vss = 1e-9*(pi*(dimensions.outerRadius-rsi)^2-Q*Su)*l; 
    
    % Magnet volume [m^3]
    Vm = 1e-9*(0.5*dimensions.rotor.alpha_m*(2*(rsi - dimensions.delta)*dimensions.rotor.dm - dimensions.rotor.dm^2)*l); 
    
    % Slot radius [mm]
    r_slot = dimensions.outerRadius - (dimensions.stator.dsy+dimensions.stator.dst/2); 
    tau_u = alpha_c * r_slot; % [mm]

    % Average length of a coil [m]
    l_c = (2*l + pi*(dimensions.stator.wst + tau_u)/2 + 2*k_ov*tau_u*(winding.y - 1))/1000; 
    
    % Total copper volume [m^3]
    Vcu = pi*dw^2/4*zQ*Q*l_c; 
    
    % Active material cost [$]
    AMC = (Vrs+Vss)*steelCost + Vcu*cuCost + Vm*magnetCost;
    
    % O_1: Active material cost [$].
    O_1 = AMC; 
        
    % O_2: Efficiency [%]
    O_2 = -designEval.efficiency; 
    
    % O_3: Torque ripple [%]
    O_3 = designEval.torque.ripple;
      
    %% Evaluate
    evaluatedChrom = [O_1, O_2, O_3];
end