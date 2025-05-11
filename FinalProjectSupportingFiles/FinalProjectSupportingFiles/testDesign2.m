clc; clear; close all;


pars = {3.8676, 1.5000, 62.1125, 22.3875, 26.7281, 25.0000};
[dm, delta, dsy, dst, wst, ast] = pars{:};

p = 2;
alpha_m = 180/p;

zQ_analytic = 317;

dso = 2;           %[mm]
dsp = 4;           %[mm]
rs = 180;          %[mm]

B_delta_hat = 0.9; % [T]
A_rms = 55e3; % [A/m]
A_hat = A_rms * sqrt(2); % [A/m]

P_rated = 50e3;
Omega_rated = 10e3 * (2*pi)/60;

ratedTorque = P_rated/Omega_rated;

Vr = ratedTorque/(B_delta_hat*A_hat);

r_rotor = 90;

dri = r_rotor - dm;

% Calculating the axial length
l = Vr/(pi*(r_rotor*1e-3)^2); % [m]
l_mm = l*1e3; % [mm]

alpha_st = ast;

kcu = 0.4;         

% Choice of wire
coil_AWG = 18;
dwire=0.324861*0.0254*exp(-0.115942*coil_AWG); % wire diameter in meters as a function of AWG
Sc = 0.25*pi*dwire^2; % [m^2]

% Peak RMS current going through the wire
J = 5; % [A_rms/mm^2]
iRMS = 5*Sc*1e6; % [A_rms]
iHat = iRMS * sqrt(2);

addpath('C:\femm42\mfiles'); % This is the FEMM path. You may need to modify this to the installation directory of FEMM on your PC.
clc;

% steps = 10:2:30;
steps = 1;

%% Stator Dimensions
dimensions.stator.dso = dso;           %[mm]
dimensions.stator.dsp = dsp;         %[mm]
dimensions.stator.dst = dst;          %[mm]
dimensions.stator.dsy = dsy;          %[mm]
dimensions.stator.alpha_st = alpha_st;     %[degree]
dimensions.stator.wst = wst;          %[mm]
dimensions.outerRadius = rs;        %[mm]
dimensions.delta = delta;              %[mm] Airgap length
dimensions.length = l_mm;             %[mm] Axial length
Q = 12;

%% Rotor Dimensions
dimensions.rotor.dri = dri;         %[mm] 
dimensions.rotor.dm=dm;            %[mm] 
dimensions.rotor.dmp=dm;            %[mm]
dimensions.rotor.alpha_m = alpha_m;     %[deg]
dimensions.rotor.OuterRadius = dimensions.outerRadius - (dimensions.stator.dsp + dimensions.stator.dst + dimensions.stator.dsy + dimensions.delta);

%% Materials 
materials.coil.name = '18 AWG';
materials.coil.conductivity = 5.77e7;   %Copper conductivity in S/m
materials.coil.AWG = 18;
materials.iron.name = 'M19_29Ga';
materials.iron.ch = 0.0186;             % Hysteresis loss coefficient for M19-29Ga in (Watt/(kg * T^2 * Hz^2)
materials.iron.ce = 6.887e-5;           % Eddy loss Current coefficientfor M19-29Ga in (Watts/(kg * T^2 * Hz)
materials.iron.sf = 0.91;               % Lamination stacking factor (nondimensional)
materials.iron.density = 7650;          % Mass density of M19 29Ga in [kg/m^3]
materials.magnet.name = 'N40';
materials.magnet.conductivity = 5.55e5; % Magnet conductivity in S/m
    
%% Define the winding structure and add parameters
winding.layers = 2; 
winding.topSlots.circuit = [{'U'}, {'W'}, {'V'}, {'U'}, {'W'}, {'V'}, {'U'}, {'W'}, {'V'}, {'U'}, {'W'}, {'V'}];
winding.topSlots.zQ = zQ_analytic.*[1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1]; % Positive number for conductors that go into the page
winding.bottomSlots.zQ = zQ_analytic.*[1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1];
winding.bottomSlots.circuit = [{'U'}, {'W'}, {'V'}, {'U'}, {'W'}, {'V'}, {'U'}, {'W'}, {'V'}, {'U'}, {'W'}, {'V'}];
winding.y = 3;
    

%% Add parameters to the settings structure
settings.RPM = 10000; % Rated Speed in RPM
settings.iHat = iHat;  % Peak current [A]
settings.phi = 0;   % Angle in [deg] (The current linkage was found in ALE15 for the same specifications, the peak is at 30 degrees, but there are two poles so 15)
settings.lowestHarmonic = 2; %Lowest Harmonic (n) for which winding factor |Kw,n| is non-zero; The max value it can take is p

%% Define the number of pole pairs
p = 2;

%% Analysis
winding_loss = NaN(size(steps));
magnet_loss = NaN(size(steps));
rotorIron_loss = NaN(size(steps));
statorIron_loss = NaN(size(steps));

for i = 1:length(steps)
    settings.steps = steps(i); % Number of steps the design must be evaluated at; Keep this an even number;

    [designEval, length] = evaluateDesign(materials, dimensions, p, winding, settings, ratedTorque);

    winding_loss(i) = designEval.loss.winding;
    magnet_loss(i) = designEval.loss.magnets;
    rotorIron_loss(i) = designEval.loss.rotorIron;
    statorIron_loss(i) = designEval.loss.statorIron;
    
    % fprintf('The avg. torque is %1.3f [Nm]; the expected value is 47.746[Nm]\n', designEval.torque.average);
    % fprintf('The torque ripple is %1.3f; the maximum value should be 0.5\n', designEval.torque.ripple);
    % fprintf('The Efficiency is %1.3f percent; the maximum value should be 90 percent\n', designEval.efficiency);
end
total_loss_1 = magnet_loss + rotorIron_loss + statorIron_loss;
total_loss_2 = total_loss_1 + winding_loss;

figure(2)
plot(steps,total_loss_1,'ro-',LineWidth=2)
xlabel('Steps')
ylabel('Total Magnet and Iron Losses (J)')
xticks(steps)
grid on
title('Plot of Total Magnet and Iron Losses vs. Number of Steps')
