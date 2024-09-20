close all
clear all

% Version for the case where a and b are independent exponentials with
% different means

% Set seed
rng(5);

% Define priors for mean and standard deviation of population
% distributions. Assuming uniform currently.

prior_mean_a_min = 0; % Upper bound on prior for mean(a).
prior_mean_a_max = 1;  % Lower bound on prior for mean(a).
prior_mean_b_min = 0;   % Upper bound on prior for mean(b).
prior_mean_b_max = 4.0;   % Lower bound on prior for mean(b).

% Define sampling parameters
nTrials = 500;         % Number of trials
bestFraction = 0.8;     % Fraction of trials to accept
bestNumber = ceil(nTrials*bestFraction); % Number of accepted trials.
samplePopDist = @sampleGamma_Gamma; % Specify the pop-level distribution


% Define synthetic data parameters
nIndiv = 10;          % Number of individuals in the dataset.
par.x0 = 1;             % Initial conditions.
par.tEnd = 1;         % Final time of model.
par.dt = 0.2;    % Time step in model.
nTime = par.tEnd/par.dt + 1;       % Number of time points.
par.sNoise = 0.1;      % Standard deviation of observation noise.
tObs = linspace(0,par.tEnd,nTime); % Vector of time points

true_mean_a = 0.1;        % True mean(a) value.
true_mean_b = 1.0;      % True mean(b) value.

true_means = [true_mean_a; true_mean_b];                        % Vector of true means.
true_theta = [true_mean_a; true_mean_b; true_mean_a^2; true_mean_b^2; 0];

% Generate sample parameters for the synthetic data.
dataParams = zeros(nIndiv,2);
dataParams = samplePopDist(true_theta, nIndiv);   % Samples from a-distribution.


% Generate synthetic data
data = evalForwardMdl(dataParams, tObs, par) + par.sNoise*normrnd(0, 1, nIndiv, nTime);


% Generate samples from prior distributions for checking
mean_a_samples = prior_mean_a_min + (prior_mean_a_max-prior_mean_a_min)*rand(nTrials,1);    % Samples of mean(a).
mean_b_samples = prior_mean_b_min + (prior_mean_b_max-prior_mean_b_min)*rand(nTrials,1);    % Samples of mean(b).

Theta_samples = [mean_a_samples, mean_b_samples, mean_a_samples.^2, mean_b_samples.^2, zeros(nTrials, 1)];

loglikelihood = zeros(nTrials,1);

% Perform trials and calculate log likelihood for each sample from the
% priors.
for iTrial = 1:nTrials
    loglikelihood(iTrial) = calcLogLik(Theta_samples(iTrial, :), data, samplePopDist, par);
end

% Sort log likelihood to find closest matches
[~,order] = sort(-loglikelihood);

% Scatter plot of estimates of mean(a) vs mean(b) for closest samples.
figure; scatter(mean_a_samples(order(1:bestNumber)),mean_b_samples(order(1:bestNumber)),15,loglikelihood(order(1:bestNumber)))
colorbar;
xline(true_means(1), 'k--')
yline(true_means(2), 'k--')
xlabel('Mean(a)'); ylabel('Mean(b)');
title('log likelihood')

