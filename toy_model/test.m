clear
close all

% Script for testing likelihood functions

% Set fixed parameter values
par.x0 = 1;
par.dt = 1;
par.sNoise = 0.01;

% Specify the pop-level distribution
samplePopDist = @sampleNormal_Gamma;

% Make up some individual parameter values for c=a+b
cIndiv = repmat( [0.4; 0.5; 0.6], 10, 1);

% Sampling time points
t = 0:5;

% Simulate synthetic data for each individual and each time point
yData = par.x0 * exp(cIndiv.*t);

% An example Theta
Theta = [0.25, 0.25, 0.01, 0.01, 0];

% Call likelihood function
tic; LLMC = calcLogLik(Theta, yData, samplePopDist, par), toc;


