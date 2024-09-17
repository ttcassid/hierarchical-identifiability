% Generate data from Friberg model for docetaxel neutrophil response
close all
clear all

%For reporducibility
rng(110);
NPatients = 5000;
% Matlab samples from Normal(mean,SD) so no need to square the SD values here. 
% normrnd(mu,sigma) generates a random number from the normal distribution with mean parameter mu and standard deviation parameter sigma.

% Define parameter distributions taken from Docetaxel response in table 4 of Friberg et al.  (2002) 
Circ0DistributionMean = 5.03; 
Circ0DistributionSD = 0.25.*sqrt( 0.15.*Circ0DistributionMean ); 
Circ0Distribution = exp( normrnd(  log(Circ0DistributionMean) ,Circ0DistributionSD,NPatients,1) ) ; 

MTTDistributionMean = 4/89.3;
MTTDistributionSD = sqrt( 0.16.*MTTDistributionMean );
kTrDistribution = exp( normrnd(  log(MTTDistributionMean) ,MTTDistributionSD,NPatients,1) ) ; 
% kTrDistribution = 4./MTTDistribution;

HSCResidenceTimeDistributionMean = 4./110.4;
HSCResidenceTimeDistributionSD = sqrt( 0.20.*HSCResidenceTimeDistributionMean );
kPfDistribution = exp( normrnd(  log(HSCResidenceTimeDistributionMean) ,HSCResidenceTimeDistributionSD,NPatients,1) ) ; 
% kPfDistribution = 1./HSCResidenceTimeDistribution;

EC50DistributionMean = 0.144; 
EC50DistributionSD = sqrt( 0.77.*EC50DistributionMean ); 
EC50Distribution = exp( normrnd(  log(EC50DistributionMean) ,EC50DistributionSD,NPatients,1) ) ; 

%No variability
GammaDistributionMean = 0.163;
GammaDistribution = GammaDistributionMean.*ones(NPatients,1);

kElDistributionMean = 1.15/24;
kElDistribution = kElDistributionMean.*ones(NPatients,1);

EMaxDistributionMean = 1;
EMaxDistribution = EMaxDistributionMean.*ones(NPatients,1);


PopulationParameters  = [Circ0Distribution, kTrDistribution, kPfDistribution, EC50Distribution,GammaDistribution,EMaxDistribution,kElDistribution];

fig = figure(1);
g1 = histogram(Circ0Distribution,'Normalization','pdf');
title('Circ0')

fig = figure(2);
g2 = histogram(kTrDistribution,'Normalization','pdf');
title('kTr')

fig = figure(3);
g3 = histogram(kPfDistribution,'Normalization','pdf');
title('kPf')

fig = figure(4);
g4 = histogram(EC50Distribution,'Normalization','pdf');
title('EC50')

% Viral dynamics parameters
fig = figure(5);
g1 = histogram( GammaDistribution,'Normalization','pdf');
title('Gamma')

fig = figure(6);
g1 = histogram( EMaxDistribution,'Normalization','pdf');
title('rhoV') 

fig = figure(6);
g1 = histogram( kElDistribution,'Normalization','pdf');
title('kEl') 
