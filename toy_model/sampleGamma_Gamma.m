function paramsIndiv = sampleGamma_Gamma(Theta, nSamples)

% Only works currently where the two individual parameters uncorrelated at
% pop level
assert(Theta(5) == 0);

% Get gamma shape and scale parameters
[sh1, sc1] = gamShapeScale(Theta(1), sqrt(Theta(3)) );
[sh2, sc2] = gamShapeScale(Theta(2), sqrt(Theta(4)) );

% Sample from the distribution (nSamples x 2 matrix where row i contains the
% two parameters for individual i)
paramsIndiv = nan(nSamples, 2);
paramsIndiv(:, 1) = gamrnd(sh1, sc1, nSamples, 1 );
paramsIndiv(:, 2) = gamrnd(sh2, sc2, nSamples, 1 );

