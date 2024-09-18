function LL = calcLogLikQuad(Theta, yData, par)

% WARNING: DON'T USE THIS FUHCTION IT ISN'T WORKING
%
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

disp('WARNING: DONT USE THE calcLogLikQuad FUNCTION IT IS NOT WORKING' )


NSD = 5;

opts = optimoptions('fmincon', 'Display', 'off');

% Extract mean and (co)variance parameters from the vector Theta
popMean = Theta(1:2);
popCov = [Theta(3), Theta(5); Theta(5), Theta(4) ];

popMin = popMean - NSD*sqrt(diag(popCov)' );
popMax = popMean + NSD*sqrt(diag(popCov)' );

[nIndiv, nPoints] = size(yData);

% Vector of time points - 3rd dimension
t = reshape( par.dt * (0:nPoints-1), 1, 1, nPoints); 


LLi = zeros(nIndiv, 1);
for iIndiv = 1:nIndiv
    yi = reshape(yData(iIndiv, :), 1, 1, nPoints);
    myLF = @(x)( -logIntegrand(x(1), x(2), t, yi, popMean, popCov, par )  );
    [~, fval] = fmincon(myLF, [Theta(1:2)], [], [], [], [], popMin, popMax, [], opts );
    LFMax = -fval;

    myF = @(a, b)( exp(logIntegrand(a, b, t, yi, popMean, popCov, par ) - LFMax)  );
    LLi(iIndiv) = LFMax + log(integral2(myF, popMin(1), popMax(1), popMin(2), popMax(2) ));
end

LL = sum(LLi);

end



function f = logIntegrand(a, b, t, yData, popMean, popCov, par)
    
    
    f1 = -1/(2*par.sNoise^2) * sum( (yData - par.x0*exp((a+b).*t)).^2, 3);

    [mRows, nCols] = size(a);
    aVec = a(:);
    bVec = b(:);

    f2 = log( reshape( mvnpdf( [aVec, bVec], popMean, popCov ), mRows, nCols) );
    f = f1 + f2;

end