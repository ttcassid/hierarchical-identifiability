close all
clear all

% Version for the case where a and b are independent exponentials with
% different means

% Set seed
rng(5);

% Define true means
true_mean_a = 0.1;        % True mean(a) value.
true_mean_b = 1.0;      % True mean(b) value.

% Define priors for mean and standard deviation of population
% distributions. Assuming uniform currently.
prior_mean_a_min = 0; % Upper bound on prior for mean(a).
prior_mean_a_max = 2.0;  % Lower bound on prior for mean(a).
prior_mean_b_min = 0;   % Upper bound on prior for mean(b).
prior_mean_b_max = 2.0;   % Lower bound on prior for mean(b).

% Define sampling parameters
nTrials = 4000;         % Number of trials
bestFraction = 0.8;     % Fraction of trials to accept
samplePopDist = @sampleExp_Exp; % Specify the pop-level distribution


% Define synthetic data parameters
nIndiv = 10;          % Number of individuals in the dataset.
par.x0 = 1;             % Initial conditions.
par.dt = 0.2;    % Time step in model.
nTime = 6;       % Number of time points.
par.sNoise = 0.1;      % Standard deviation of observation noise.



tObs = par.dt*(0:nTime-1); % Vector of observation time points


true_means = [true_mean_a; true_mean_b];                        % Vector of true means.
true_theta = [true_mean_a; true_mean_b; true_mean_a^2; true_mean_b^2; 0];       % True 'theta' vector

% Generate sample parameters for the synthetic data.
dataParams = samplePopDist(true_theta, nIndiv);   % Samples from a-distribution.


% Generate synthetic data
logData = evalForwardMdl(dataParams, tObs, par) + par.sNoise*normrnd(0, 1, nIndiv, nTime);


% Generate samples from prior distributions for checking
mean_a_samples = prior_mean_a_min + (prior_mean_a_max-prior_mean_a_min)*rand(nTrials,1);    % Samples of mean(a).
mean_b_samples = prior_mean_b_min + (prior_mean_b_max-prior_mean_b_min)*rand(nTrials,1);    % Samples of mean(b).
Theta_samples = [mean_a_samples, mean_b_samples, mean_a_samples.^2, mean_b_samples.^2, zeros(nTrials, 1)];

loglikelihood = zeros(nTrials,1);

% Perform trials and calculate log likelihood for each sample from the
% priors.
parfor iTrial = 1:nTrials
    loglikelihood(iTrial) = calcLogLik(Theta_samples(iTrial, :), logData, samplePopDist, par);
end

% Sort log likelihood to find closest matches
[~,order] = sort(-loglikelihood);

bestNumber = ceil(nTrials*bestFraction); % Number of accepted trials.

% Scatter plot of estimates of mean(a) vs mean(b) for closest samples.
figure; scatter(mean_a_samples(order(1:bestNumber)),mean_b_samples(order(1:bestNumber)),15,loglikelihood(order(1:bestNumber)))
colorbar;
xline(true_means(1), 'k--')
xline(true_means(2), 'k--')
yline(true_means(1), 'k--')
yline(true_means(2), 'k--')

xlabel('Mean(a)'); ylabel('Mean(b)');
title('log likelihood')

