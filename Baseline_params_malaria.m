%% Malaria parameters taken from Qu and Patterson et al 2023
P.bm = 0.6;
P.bh = 5;

P.betaM = 0.25;
P.betaD = 0.35;
P.betaA = 0.03;

P.rD = 1/180;
P.rA = 1/360;

P.h = 1/15;
P.muD = 0;

P.cS = 0.75;
P.cE = 0.1;
P.cA = 0.1;
P.cD = 0.05;

P.de = 5*365;
P.gamma = 10;

% Fixed immunity parameters
P.psi = 0.5; % descreasing function in immunity
P.rho = 0.5; % descreasing 
P.phi = 0.5; % increasing 

% Simple birth/death parameters (lifespan = 80)
P.gH = 1;
P.muH = 1/(80*365);
