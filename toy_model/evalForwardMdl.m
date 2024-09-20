function xt = evalForwardMdl(paramsIndiv, tObs, par )

cIndiv = sum(paramsIndiv, 2);                                    % a+b

xt = par.x0*exp(cIndiv.*tObs);

