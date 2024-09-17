% Virutal population simulation to generate neutrophil response data
% to 63 days of Zalypsis (with parameters taken from Craig et al. (2016),
% Cassidy et al. 2021). 
% Implemented the Friberg model following their parameterization (including
% the misspecification of transit rates (see de Souza et al. JPKPD) 
clear all
close all
%% Choose dose size and load virtual population
FribergVPopStruct = load('VirtualPopulationParameters_FribergModel_17Sept24'); %   load('TwoStrainMutation_CombobNAB_ViralDynamics_VRC07Param_23Dec2022.mat'); %
FribergHBVVPop = cell2mat(struct2cell(FribergVPopStruct));
%% PK Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters for the chemotherapy (Zalypsis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PA.Cl=43.7*24; %Rate of clearance
PA.V1=32.7; %Volume in central compartment
PA.Q2=123*24; %Amount of transit (second compartment)
PA.V2=162; %Volume of second compartment
PA.Q3=11.3*24; %Amount of transit (third compartment)
PA.V3=388; %Volume of third compartment
PA.Q4=62.3*24;%Amount of transit (fourth compartment)
PA.V4=23.9;%Volume of fourth compartment
PA.k10=PA.Cl/PA.V1; %Rate of elimination
PA.k12=PA.Q2/PA.V1;%Rate of transit from compartment 1 to 2
PA.k21=PA.Q2/PA.V2;%Rate of transit from compartment 2 to 1
PA.k13=PA.Q3/PA.V1;%Rate of transit from compartment 1 to 3
PA.k31=PA.Q3/PA.V3;%Rate of transit from compartment 3 to 1
PA.k24=PA.Q4/PA.V2;%Rate of transit from compartment 2 to 4
PA.k42=PA.Q4/PA.V4;%Rate of transit from compartment 4 to 2
PA.BSA=1.723;%Average body surface area


%% Set up Bootstrapping for trial size
NPatients =  10 ; % length(FribergHBVVPop(:,1));
NTrials = 1;
%% Set up virtual clinical trial calculate proportion of virus
% Simulation time
t0 = 0;
tf =  65;
TotalTime = [t0 tf];

%% Chemo dosing
PA.CycleLength =  21; %Set the treatment schedule length in days
PA.TreatmentDuration = 1; %Number of doses per cycle
PA.TreatmentPerCycle = 1;  %Treatments per cycle, including the initial dose
PA.AdminNumber = ceil(tf./PA.CycleLength) ;
PA.TimeAdmin = 1/24; %  %Administration time of treatment
PA.StartTime = t0; %100; %  %start of treatment .

PA.DoseChemo = 4000; % Chemotherapy dose (mcg/m^2) %Administer a dose of 0 chemo
PA.TotalDose =PA.DoseChemo*PA.BSA;%Chemotherapy dose (mcg)
PA.Admin = ones(1,PA.AdminNumber).*PA.TotalDose;

PA.TreatmentStartTime = PA.StartTime + [0:PA.AdminNumber-1].*PA.CycleLength ; %Days when treatment occurs
PA.TreatmentEndTime = PA.TreatmentStartTime+PA.TimeAdmin;

%% Sampling times for neutrophil concentrations:
SampleTimes = [t0:3:tf] ;
NSamples = length(SampleTimes); 
%% Set up results for each virtual population
VirtualPatientNeutrophilConcentration = zeros(NSamples,NPatients);
%% Set up virtual clinical trial calculate proportion of virus
tic

for ii = 1: NPatients  % For each participant in the bootstrapped trial
    
    %% Load individual virtual patient parameter estimates
    PatientNumberParameter = ii; % SampleVpatientIndex(jj,ii);
    X = FribergHBVVPop(PatientNumberParameter,1:end);
    %% Update model parameters
    % PopulationParameters  = [Circ0Distribution, kTrDistribution, kPfDistribution, EC50Distribution,GammaDistribution,EMaxDistribution,kElDistribution];

    PA.circ0 = X(1);
    PA.ktr = X(2);
    PA.kprol = X(3);
    PA.EC50 =  X(4);
    PA.gamma = X(5);
    PA.Emax = X(6);;
    PA.kcirc = X(7);
    PA.N = 3; 
    
    % Parameters Erlang
    PA.tau = PA.N ./ (PA.ktr);
    
    %% Initial condtions 
    CircIC = PA.circ0;
    ProlifIC = PA.kcirc ./ PA.kprol .* PA.circ0;
    ICTransitCompartments = PA.kcirc ./  PA.ktr .* PA.circ0 .* ones(1,PA.N);
    Chemo1IC = 0; % PA.TotalDose;
    Chemo2IC = 0;
    Chemo3IC = 0;
    Chemo4IC = 0;
    IC = [ProlifIC,ICTransitCompartments,CircIC,Chemo1IC,Chemo2IC,Chemo3IC,Chemo4IC];
    %% Solve the ODE systems 
    
    [sol_myelosuppression_model] = Myelosuppression_Model(TotalTime, IC, PA); % Solve until end of experiment 
    NeutrophilSamples = deval(sol_myelosuppression_model, SampleTimes,PA.N+2);
    VirtualPatientNeutrophilConcentration(:,ii) = NeutrophilSamples; 
end
toc
VirtualPatientNeutrophilConcentration = [SampleTimes',VirtualPatientNeutrophilConcentration];
% quantile(VirtualPatientNeutrophilConcentration,[0.05,0.10,0.50,0.90,0.95]);

% quantile(MeanRNAChange,[0.05,0.10,0.50,0.90,0.95])

t = datetime('now','Format','dMMMyy');
S = char(t);
filename = [S,'_FullVPopCohortSize_FribergDoseResponse'];
save(filename)


function Dose = ChemoDose(PA,t); %IV administration of chemotherapy
DoseVec = zeros(1,PA.AdminNumber); %Create a vector of doses for each administration
for nn = 1:PA.AdminNumber
    if PA.TreatmentStartTime(nn)  < t && t < PA.TreatmentEndTime(nn)
        DoseVec(nn) = (PA.Admin(nn))./(PA.TimeAdmin);
    else
        DoseVec(nn) = 0;
    end
end
Dose = ones(1,PA.AdminNumber)*DoseVec';
end

function [sol] = Myelosuppression_Model(TotalTime, IC, PA)
opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',1e-2);
sol = ode45(@Myelosuppression_Dynamics,TotalTime,IC,opts);
    function dydt = Myelosuppression_Dynamics(t, y);
        dydt(1) = PA.kprol * y(1) * (1-PA.Emax.*y((PA.N+2)+1) ./ (PA.EC50+y((PA.N+2)+1) )) * (PA.circ0 ./ y(PA.N+2))^PA.gamma - PA.kprol .* y(1);
        dydt(2) = PA.kprol*y(1) - PA.ktr*y(2);
        for ii = 3:PA.N+1
            dydt(ii) = PA.ktr*(y(ii-1) - y(ii));
        end
        dydt(PA.N+2) = PA.ktr * y(PA.N+2-1) - PA.kcirc * y(PA.N+2);
        dydt((PA.N+2)+1) = ChemoDose(PA,t) - (PA.k10+PA.k12+PA.k13)*y((PA.N+2)+1) + PA.k21*y((PA.N+2)+2 ) + PA.k31*y((PA.N+2)+3 ); %Chemo 1 compartment DE
        dydt((PA.N+2)+2) = -(PA.k21+PA.k24)*y((PA.N+2)+2)+ PA.k12*y((PA.N+2)+1) + PA.k42*y((PA.N+2)+4) ; %Chemo 2 compartment
        dydt((PA.N+2)+3) = -PA.k31*y((PA.N+2)+3)+PA.k13*y((PA.N+2)+1); %Chemo 3 compartment
        dydt((PA.N+2)+4) = -PA.k42*y((PA.N+2)+4)+PA.k24*y((PA.N+2)+2);%Chemo 4 compartment 
        dydt = dydt';
    end
end











