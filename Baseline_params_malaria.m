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
P.psi = 0.5;
P.rho = 0.5;
P.phi = 0.5;

% Simple birth/death parameters (lifespan = 80)
P.gH = 1;
P.muH = 1/(80*365);

%{
%% Wolbachia parameters taken from Qu, Wu, Hyman 2022
P.vw = 1;
P.vu = 1-P.vw;
P.bf = 0.5;
P.Kf = 3e5;
P.ci = 1; % new parameter

P.phiU = 7.0;
P.phiW = 5.7;
P.mufu = 1/26.25;
P.mufw = 1/24.25;

%% New parameters linking malaria and Wolbachia
P.sigma = 1/10;
P.alpha = 0;

%}