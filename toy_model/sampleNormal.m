function paramsIndiv = sampleNormal(Theta, nSamples)

% Extract mean and (co)variance parameters from the vector Theta
popMean = Theta(1:2);
popCov = [Theta(3), Theta(5); Theta(5), Theta(4) ];

% Sample from the distribution (nSamples x 2 matrix where row i contains the
% two parameters for individual i)
paramsIndiv = mvnrnd(popMean, popCov, nSamples );