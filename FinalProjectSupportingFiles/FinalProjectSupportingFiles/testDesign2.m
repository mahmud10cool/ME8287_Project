%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a test script for HW 3 Part 2). This is an incomplete script and
% you need to complete it as per the instructions provided.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
clear all;
close all;
addpath('C:\femm42\mfiles'); % This is the FEMM path. You may need to modify this to the installation directory of FEMM on your PC.
clc;

%% Stator Dimensions
dimensions.stator.dso = 5;           %[mm]
dimensions.stator.dsp = 7.5;         %[mm]
dimensions.stator.dst = 57;          %[mm]
dimensions.stator.dsy = 35;          %[mm]
dimensions.stator.alpha_st = 55;     %[degree]
dimensions.stator.wst = 50;          %[mm]
dimensions.outerRadius = 200;        %[mm]
dimensions.delta = 3.5;              %[mm] Airgap length
dimensions.length = 600;             %[mm] Axial length
Q = 6;

%% Rotor Dimensions
dimensions.rotor.dri = 50;         %[mm] 
dimensions.rotor.dm=8;            %[mm] 
dimensions.rotor.dmp=6;            %[mm]
dimensions.rotor.alpha_m = 80;     %[deg]
dimensions.rotor.OuterRadius = dimensions.outerRadius - (dimensions.stator.dsp + dimensions.stator.dst + dimensions.stator.dsy + dimensions.delta);

%% Materials 
materials.coil.name = '18 AWG';
materials.coil.conductivity = 5.77e7; %Copper conductivity in S/m
materials.coil.AWG = 18;
materials.iron.name = 'M19_29Ga';
materials.iron.ch = 0.0186;  % Hysteresis loss coefficient for M19-29Ga in (Watt/(kg * T^2 * Hz^2)
materials.iron.ce = 6.8874e-5; % Eddy loss Current coefficient for M19-29Ga in (Watts/(kg * T^2 * Hz)
materials.iron.sf = 0.91;  % Lamination stacking factor (nondimensional)
materials.iron.density = 7650; %Mass density of M19 29Ga in [kg/m^3]
materials.magnet.name = 'Recoma35E';
materials.magnet.conductivity = 1.11e6; %Magnet conductivity in S/m

%% Define the winding structure and add parameters
dwire=0.324861*0.0254*exp(-0.115942*materials.coil.AWG); % wire diameter in meters as a function of AWG
Acond = 0.25*pi*dwire^2;
iHat = 5*1e6*Acond*sqrt(2);
A_hat = 25e3;
N = Q/3; % Because it is a double layer winding
r_airgap = 1e-3*(dimensions.rotor.OuterRadius);

zQ = round((pi*A_hat*r_airgap)/(3*N*0.866*iHat)); % This is 261 turns
% zQ = 252; % This is just to match the given torque values.
winding.layers = 2; 
winding.topSlots.circuit = [{'U'}, {'V'}, {'W'},{'U'}, {'V'}, {'W'}];
winding.topSlots.zQ = zQ.*[1, 1, 1, 1, 1, 1]; % Positive number for conductors that go into the page
winding.bottomSlots.zQ = zQ.*[-1, -1, -1, -1, -1, -1];
winding.bottomSlots.circuit = [{'W'}, {'U'}, {'V'},{'W'}, {'U'}, {'V'}];
winding.y = 1;

%% Add parameters to the settings structure
settings.RPM = 3600; % Rated Speed in RPM
settings.steps = 12; % Number of steps the design must be evaluated at; Keep this an even number;
settings.iHat = iHat;  % Peak current [A]
settings.phi = 30;   % Angle in [deg] (The current linkage was found in ALE15 for the same specifications, the peak is at 30 degrees, but there are two poles so 15)
settings.lowestHarmonic = 2; %Lowest Harmonic (n) for which winding factor |Kw,n| is non-zero; The max value it can take is p

%% Define the number of pole pairs
p = 2;

%% Analysis
[designEval] = evaluateDesign(materials, dimensions, p, winding, settings);

fprintf('The avg. torque is %1.3f [Nm]; the expected value is 398.058[Nm]\n', designEval.torque.average);
fprintf('The torque ripple is %1.3f; the expected value is 0.130\n', designEval.torque.ripple);
fprintf('The Efficiency is %1.3f percent; the expected value is 86.635 percent\n', designEval.efficiency);

