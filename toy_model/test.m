clear
close all

% Script for testing likelihood functions

% Set fixed parameter values
par.x0 = 1;
par.dt = 0.2;
par.sNoise = 0.01;
nPoints = 6;

% Specify the pop-level distribution
samplePopDist = @sampleNormal_Gamma;

% Make up some individual parameter values for c=a+b
cIndiv = [   1.6086;     1.7271;    1.1718;    2.0941;    1.5350;    1.5476;    1.4874;    1.1833;    1.3962;    0.6735];

% Sampling time points
t = par.dt*(0:nPoints-1);

% Simulate synthetic data for each individual and each time point
yData = par.x0 * exp(cIndiv.*t);

% An example Theta
Theta = [0, 1.5, 0.01, 0.01, 0];

% Call likelihood function
tic; LLMC = calcLogLik(Theta, yData, samplePopDist, par), toc;


