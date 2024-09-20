function logxt = evalForwardMdl(paramsIndiv, tObs, par )

cIndiv = sum(paramsIndiv, 2);                                    % a+b

logxt = log(par.x0) + cIndiv.*tObs;

