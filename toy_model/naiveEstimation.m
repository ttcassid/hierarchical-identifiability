close all
clear all

% Set seed
rng(5);

% Define priors for mean and standard deviation of population
% distributions. Assuming uniform currently.

prior_mean_a_min = -1; % Upper bound on prior for mean(a).
prior_mean_a_max = 1;  % Lower bound on prior for mean(a).
prior_sd_a_min = 0;     % Upper bound on prior for std(a).
prior_sd_a_max = 1;     % Lower bound on prior for std(a).
prior_mean_b_min = 0;   % Upper bound on prior for mean(b).
prior_mean_b_max = 3;   % Lower bound on prior for mean(b).
prior_sd_b_min = 0;     % Upper bound on prior for std(b).
prior_sd_b_max = 3;     % Lower bound on prior for std(b).
prior_cov_min = 0;      % Upper bound on prior for cov(a,b).
prior_cov_max = 0;      % Lower bound on prior for cov(a,b).

% Define sampling parameters
nTrials = 500;         % Number of trials
bestFraction = 0.1;     % Fraction of trials to accept
bestNumber = ceil(nTrials*bestFraction); % Number of accepted trials.
samplePopDist = @sampleNormal_Gamma; % Specify the pop-level distribution


% Define synthetic data parameters
nIndiv = 10;          % Number of individuals in the dataset.
par.x0 = 1;             % Initial conditions.
par.tEnd = 0.1;         % Final time of model.
par.dt = par.tEnd/4;    % Time step in model.
nTime = 1/par.dt;       % Number of time points.
par.sNoise = 0.01;      % Standard deviation of observation noise.
t = linspace(0,par.tEnd-par.dt,nTime); % Vector of time points

true_mean_a = 0;        % True mean(a) value.
true_sd_a = 0.3;        % True std(a) value.
true_mean_b = 1.5;      % True mean(b) value.
true_sd_b = 2;          % True std(b) value.
true_cov = 0;           % True cov(a,b) value.

true_means = [true_mean_a; true_mean_b];                        % Vector of true means.
true_sds = [true_sd_a; true_sd_b];                              % Vector of true stds.
true_covMatrix = [true_sd_a, true_cov; true_cov, true_sd_b];    % True covariance matrix.

% Generate sample parameters for the synthetic data.
dataParams = zeros(nIndiv,2);
dataParams(:,1) = normrnd(true_means(1),true_sds(1),nIndiv,1);  % Samples from a-distribution.
dataParams(:,2) = gamrnd(true_means(2),true_sds(2),nIndiv,1);   % Samples from b-distribution.
cIndiv = sum(dataParams, 2);                                    % a+b

% Generate synthetic data
data = par.x0*exp(cIndiv*t);

% Generate samples from prior distributions for checking
mean_a_samples = prior_mean_a_min + (prior_mean_a_max-prior_mean_a_min)*rand(nTrials,1);    % Samples of mean(a).
sd_a_samples = prior_sd_a_min + (prior_sd_a_max-prior_sd_a_min)*rand(nTrials,1);            % Samples of std(a).
mean_b_samples = prior_mean_b_min + (prior_mean_b_max-prior_mean_b_min)*rand(nTrials,1);    % Samples of mean(b).
sd_b_samples = prior_sd_b_min + (prior_sd_b_max-prior_sd_b_min)*rand(nTrials,1);            % Samples of std(b).
cov_samples = prior_cov_min + (prior_cov_max-prior_cov_min)*rand(nTrials,1);                % Samples of cov(a,b).

loglikelihood = zeros(nTrials,1);

% Perform trials and calculate log likelihood for each sample from the
% priors.
for iTrial = 1:nTrials
    iTrial
    Theta = [mean_a_samples(iTrial);mean_b_samples(iTrial);sd_a_samples(iTrial);sd_b_samples(iTrial);cov_samples(iTrial)];
    loglikelihood(iTrial) = calcLogLik(Theta, data, samplePopDist, par);
end

% Sort log likelihood to find closest matches
[~,order] = sort(-loglikelihood);

% Scatter plot of estimates of mean(a) vs mean(b) for closest samples.
figure; scatter(mean_a_samples(order(1:bestNumber)),mean_b_samples(order(1:bestNumber)),15,loglikelihood(order(1:bestNumber)))
xlabel('Mean(a)'); ylabel('Mean(b)');