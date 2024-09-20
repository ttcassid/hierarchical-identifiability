function paramsIndiv = sampleExp_Exp(Theta, nSamples)


% Sample from the distribution (nSamples x 2 matrix where row i contains the
% two parameters for individual i)
paramsIndiv = exprnd( repmat([Theta(1), Theta(2)], nSamples, 1)); 

