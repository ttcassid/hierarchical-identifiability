function LL = calcLogLik(Theta, yData, par)

% Function to calculate the log likelihood for population-level parameters
% in the vector Theta
% Theta = [ mean(a), mean(b), var(a), var(b), cov(a,b) ]
% 
% USAGE: LL = calcLogLik(Theta, yData, par)
% INPUTS: Theta (as defined above)
%         yData - matrix of data whose (i,j,) element is the observation or individual i and time point j
%         par - structure of fixed parameters with following fields
%         par.x0 = initial condition  for x(0)
%         par.dt = time step between observations
%         par.sNoise - SD of observation noise
 


MCMC_SAMPLES = 1e5;

% Extract mean and (co)variance parameters from the vector Theta
popMean = Theta(1:2);
popCov = [Theta(3), Theta(5); Theta(5), Theta(4) ];

[nIndiv, nPoints] = size(yData);

% Total number of required samples from pop-level dist required for Monte Carlo intergation
nSamples = nIndiv * MCMC_SAMPLES;

paramsIndiv = mvnrnd(popMean, popCov, nSamples );

% All that matters for the likelihood function is c = a+b
cIndiv = sum(paramsIndiv, 2);

% Reshape into a matrix whose (i,j) element is the jth MCMC sample for the ith individual 
% ith individual
cIndiv = reshape(cIndiv, nIndiv, MCMC_SAMPLES);

% Vector of time points (in 3rd dimension for summing over)
t = reshape( par.dt * (0:nPoints-1), 1, 1, nPoints ); 

% Reshape data matrix so that time points are in the 3rd dimension to match
% t
yReshaped = reshape(yData, nIndiv, 1, nPoints);

% Calculate (log) probabilities and sum over time points 
Pij = sum((  yReshaped  - par.x0*exp(cIndiv.* t) ).^2, 3);

LLi = logsumexp(-1/par.sNoise^2 * Pij, 2);

LL = sum(LLi);





