clear;
close all;
clc;

%% Define variables
numVariables = 6; % Number of free variables

% % variables: [dm, delta, dsy, dst, wst, ast] 
boundsLo = [0.25, 1, 5,5,5,1.8]; % Lower bounds [mm,mm,mm,mm,mm,deg]
boundsHi = [5,5,50,100,100,28];% Upper Bounds [mm,mm,mm,mm,mm,deg]

%% Define the constraint/s
A =[]; % Linear inequality constraints 
b =[]; % Linear inequality constraints 
Aeq = []; % Linear equality constraints 
beq = []; % Linear equality constraints 

%% Set the number of objectives and the function that evaluates the objectives
numObj =3; %Number of objectives
evalObjectives = @evaluateObjectives;% This will be the MATLAB function that evaluates each generation
constraint = @evaluateConstraints;% This will be the MATLAB function that evaluates constraints

%% Set the optimization options
generations = 5; % Set the number of generations
population = 50; % Set the population size (This is the number of individuals per generation)

%% Run the optimization
if isfile('population1.mat')
    disp('Starting from the saved population');
    load('population1.mat', 'savedPopulation');
    initialPopulation = savedPopulation;
else
    disp('Starting from the analytic design from Q5:');
    initialPopulation = [3.8676, 1.5000, 45.5310, 38.9690, 26.7281, 25.0000];
end

options = optimoptions(@gamultiobj, ...
    'PlotFcn', @gaplotpareto, ...
    'MaxGenerations', generations, ...
    'PopulationSize', population, ...
    'UseParallel', true, ...
    'InitialPopulationMatrix', initialPopulation, ...
    'OutputFcn', @displayGeneration);

% Run the optimization
[x,Fval,exitFlag,Output,population,scores] = gamultiobj(evalObjectives,numVariables,A, ...
    b,Aeq,beq,boundsLo,boundsHi,constraint,options);

% Final population to resume later
savedPopulation = population;
save('population2.mat', 'savedPopulation');

% Scores
save('scores_gen2.mat', 'scores');

% Pareto front
save('par_front_gen2.mat', 'Fval');


