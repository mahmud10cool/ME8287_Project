%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a test script for HW 3 Part 2). This defines the materials,
% dimensions, excitation and windings and calls the relevant functions
% to evaluate a candidate design.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
clear all;
close all;
addpath('C:\femm42\mfiles'); % This is the FEMM path. You may need to modify this to the installation directory of FEMM on your PC.
clc;

%% Stator Dimensions
dimensions.stator.dso = 2;           %[mm]
dimensions.stator.dsp = 4.5;         %[mm]
dimensions.stator.dst = 60;          %[mm]
dimensions.stator.dsy = 35;          %[mm]
dimensions.stator.alpha_st = 25;     %[degree]
dimensions.stator.wst = 25;          %[mm]
dimensions.outerRadius = 200;        %[mm]
dimensions.delta = 3.5;              %[mm] Airgap length
dimensions.length = 600;             %[mm] Axial length
Q=6;

%% Rotor Dimensions
dimensions.rotor.dri = 40;         %[mm] 
dimensions.rotor.dm=12;            %[mm] 
dimensions.rotor.dmp=5;            %[mm]
dimensions.rotor.alpha_m = 60;     %[deg]
dimensions.rotor.OuterRadius = dimensions.outerRadius - (dimensions.stator.dsp + dimensions.stator.dst + dimensions.stator.dsy + dimensions.delta);
p = 2;                             % pole pairs

%% Materials and Windings
materials.coil.name = '18 AWG';
materials.coil.conductivity = 5.77e7; %Copper conductivity in S/m
materials.coil.AWG = 18;
materials.iron.name = 'M19_29Ga';
materials.iron.ch = 0.0186;  % Hysteresis loss coefficient for M19-29Ga in (Watt/(kg * T^2 * Hz^2)
materials.iron.ce = 6.8874e-5;   % Eddy Current loss coefficient for M19-29Ga in (Watts/(kg * T^2 * Hz)
materials.iron.sf = 0.91;  % Lamination stacking factor (nondimensional)
materials.iron.density = 7650; %Mass density of M19 29Ga in [kg/m^3]
materials.magnet.name = 'Recoma35E';
materials.magnet.conductivity = 1.11e6; %Magnet conductivity in S/m

%% Winding
zQ = 80;                    
winding.layers = 2; 
winding.topSlots.circuit = [{'U'}, {'W'}, {'V'},{'U'}, {'W'}, {'V'},{'U'}, {'W'}, {'V'},{'U'}, {'W'}, {'V'}];
winding.topSlots.zQ = zQ.*[1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1]; % Positive number for conductors that go into the page
winding.bottomSlots.zQ = zQ.*[-1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1];
winding.bottomSlots.circuit = [{'W'}, {'V'}, {'U'},{'W'}, {'V'}, {'U'},{'W'}, {'V'}, {'U'},{'W'}, {'V'}, {'U'}];
winding.y = 2;

%% Settings
settings.RPM = 3600; % Rated Speed in RPM
settings.phi = 30;   % Angle in [deg]
settings.iHat = 18;  % Peak current [A]
settings.lowestHarmonic = 2; %Lowest Harmonic (n) for which winding factor |Kw,n| is non-zero; The max value it can take is p
settings.steps = 18; % Number of steps the design must be evaluated at; Let this be an even number

%% Analysis
[designEval] = evaluateDesign(materials, dimensions, p, winding, settings);

%% Results
fprintf('The avg. torque is %1.3f [Nm]; the expected value is 722.154[Nm]\n', designEval.torque.average);
fprintf('The torque ripple is %1.3f; the expected value is 0.455\n', designEval.torque.ripple);
fprintf('The Efficiency is %1.3f percent; the expected value is 94.449 percent\n', designEval.efficiency);
