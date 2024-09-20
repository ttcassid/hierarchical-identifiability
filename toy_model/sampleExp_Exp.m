function paramsIndiv = sampleExp_Exp(Theta, nSamples)


% Sample from the distribution (nSamples x 2 matrix where row i contains the
% two parameters for individual i)
paramsIndiv = exprnd( repmat(Theta, nSamples, 1)); 

