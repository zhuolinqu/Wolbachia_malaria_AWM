%% Malaria parameters taken from Qu and Patterson et al 2023 (except as noted)
P.bm = 0.6;
P.bh = 5;

P.betaM = 0.25;
P.betaD = 0.35;
P.betaA = 0.03;

P.rD = 1/33.5; % Updated based on treatment assumption using malaria therapy data % 1/180; original Qu and Patterson et al
P.rA = 1/85; % Updated  using malaria therapy data  % 1/360; original Qu and Patterson et al

P.h = 1/26; % Updated based on malaria therapy data % 1/15; original Qu and Patterson et al
P.muD = 0;

P.cS = 0.75;
P.cE = 0.1;
P.cA = 0.1;
P.cD = 0.05;

P.de = 5*365;
P.gamma = 10;

% Fixed immunity parameters
P.psi = 0.5; % increasing function in immunity
P.rho = 0.5; % increasing 
P.phi = 0.5; % increasing 

% Simple birth/death parameters (lifespan = 80)
P.gH = 1;
P.muH = 1/(80*365);
