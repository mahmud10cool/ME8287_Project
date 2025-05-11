clear;
close all;
clc;

%% Define variables
numVariables = 6; % Number of free variables

% creat first generation close to the analytic design
% variables: [dm, delta, dsy, dst, wst, ast] 
% boundsLo = [4.2,1.9,48.9,49.1,21.3,23]; % Lower bounds [mm,mm,mm,mm,mm,deg]
% boundsHi = [4.4,2.1,49.99,51.1,23.3,25];% Upper Bounds [mm,mm,mm,mm,mm,deg]

% % variables: [dm, delta, dsy, dst, wst, ast] 
boundsLo = [0.25, 1, 5,5,5,1.8]; % Lower bounds [mm,mm,mm,mm,mm,deg]
boundsHi = [5,5,50,100,100,28];% Upper Bounds [mm,mm,mm,mm,mm,deg]

%% Define the constraint/s
A =[]; % Linear inequality constraints (We are not using it for this problem)
b =[]; % Linear inequality constraints (We are not using it for this problem)
Aeq = []; % Linear equality constraints (We are not using it for this problem)
beq = []; % Linear equality constraints (We are not using it for this problem)

%% Set the number of objectives and the function that evaluates the objectives
numObj =3; %Number of objectives
evalObjectives = @evaluateObjectives;% This will be the MATLAB function that evaluates each generation
constraint = @evaluateConstraints;% This will be the MATLAB function that evaluates constraints

%% Set the optimization options
generations = 50; % Set the number of generations
population = 50; % Set the population size (This is the number of individuals per generation)
% options = optimoptions(@gamultiobj,'PlotFcn',@gaplotpareto, 'MaxGenerations',generations, 'PopulationSize',...
%                        population,'UseParallel',true);

%% Run the optimization
% [x,Fval,exitFlag,Output,population,scores] = gamultiobj(evalObjectives,numVariables,A, ...
%     b,Aeq,beq,boundsLo,boundsHi,constraint,options);
% save("HW4SimulationResults.mat")

if isfile('finalPop.mat')
    disp('Resuming from saved population...');
    load('finalPop.mat', 'savedPopulation');
    initialPopulation = savedPopulation;
else
    disp('Starting from the analytic design');
    initialPopulation = [3.8676, 1.5, 62.1125, 22.3875, 26.7281, 25.0000];
end

options = optimoptions(@gamultiobj, ...
    'PlotFcn', @gaplotpareto, ...
    'MaxGenerations', generations, ...
    'PopulationSize', population, ...
    'UseParallel', true, ...
    'InitialPopulationMatrix', initialPopulation);

% Run the optimization
[x,Fval,exitFlag,Output,population,scores] = gamultiobj(evalObjectives,numVariables,A, ...
    b,Aeq,beq,boundsLo,boundsHi,constraint,options);

% save final population to resume later
savedPopulation = population;
save('finalPop.mat', 'savedPopulation');

% save scores
save('scores_gen1.mat', 'scores');

% save pareto front
save('pareto_front_gen1.mat', 'Fval');


