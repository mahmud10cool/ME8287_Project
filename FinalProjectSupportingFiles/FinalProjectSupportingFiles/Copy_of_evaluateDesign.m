%% Homework 4
% Date: 04/24/2025

function evaluatedChrom = evaluateDesign(parameter)
% EVALUATEDESIGN evaluates a design (candidate design) and computes
% the objectives. Any constraints/constants if needed, can be defined local to this
% function. FEA runs if any, will be performed here.

%% Initialization
% Get the dimensions from the chromosome.
B_delta = parameter(1);
A_hat = parameter(2);
delta = parameter(3);
r = parameter(4);
Omega = parameter(5);
   
%% Compute Objectives.
P_rated = 100e3;

T_rated = P_rated/(Omega);
l = T_rated/(B_delta*A_hat*pi*r^2);

Vr = pi*r^2*l;

% T = Vr*B_delta*A_hat;

% Power_out = T*Omega;

% Objectives
% Minimize losses
rho = 7850;
Vs = pi*((1.5*r)^2-r^2)*l;
Ce = 6.88e-5;
Chy = 0.0186;
C_Omega = 0.002;
f = (Omega)/(2*pi);
P_loss = rho*Vs*(Ce*f^2*B_delta^2 + Chy*f*B_delta^2) + C_Omega*Vs*A_hat^2;

% Another way to optimize would be to maximize efficiency
% Efficiency = 100*(Power_out)/(Power_out+P_loss);

% Minimize torque ripple
T_ripple = 1/(20*(r*delta)^0.25)+0.25*sin(4000*pi*delta);

% Minimize active material cost
Pr = 1000*(9.85+5000*delta);
Ps = 6750;
C = Pr*Vr + Ps*Vs;
  
%% Return the Power Loss, Torque Ripple and Cost
evaluatedChrom = [C, T_ripple, P_loss];

% This is the maximizing the efficiency way
% evaluatedChrom = [C, T_ripple, -Efficiency];

end