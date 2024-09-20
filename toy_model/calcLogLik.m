function LL = calcLogLik(Theta, logData, samplePopDist, par)

% Function to calculate the log likelihood for population-level parameters
% in the vector Theta
% Theta = [ mean(a), mean(b), var(a), var(b), cov(a,b) ]
% 
% USAGE: LL = calcLogLik(Theta, yData, par)
% INPUTS: Theta (as defined above)
%         logData - matrix of data whose (i,j,) element is the logged observation or individual i and time point j
%         samplePopDist - a function to return samples of individual parameters from the population lervel distribution
%                       - called with input arguments Theta and nSamples
%         par - structure of fixed parameters with following fields
%         par.x0 = initial condition  for x(0)
%         par.dt = time step between observations
%         par.sNoise - SD of observation noise
 


MCMC_SAMPLES = 1e5;

[nIndiv, nPoints] = size(logData);


% Sample individual parameters from the pop level distribution
paramsIndiv = samplePopDist(Theta, MCMC_SAMPLES);

% Vector of time points
tObs = par.dt * (0:nPoints-1); 

% Evaluate forward model under individual parameters
logxt = evalForwardMdl(paramsIndiv, tObs, par);

% Reshape so time is on the 3rd dimension for summing over
logxt = reshape(logxt, 1, MCMC_SAMPLES, nPoints );


% Reshape data matrix so that time points are in the 3rd dimension to match
% t
logDataReshaped = reshape(logData, nIndiv, 1, nPoints);

% Calculate (log) probabilities and sum over time points 
Pij = sum(( logDataReshaped - logxt).^2, 3);

LLi = logsumexp(-1/par.sNoise^2 * Pij, 2);

LL = sum(LLi);





