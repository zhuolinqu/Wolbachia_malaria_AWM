%% Malaria parameters taken from Qu and Patterson et al 2023 (except as noted)
P.bm = 0.6;
P.bh = 5;

P.betaM = 0.2;
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

% Simple birth/death parameters (ave lifespan = 80)
P.gH = 2.5;
P.muH = 1/(60*365);

% sigmoid parameters, fitted, taken from malaria simulation project
x = [2.567957971786876   2.487540758554113   3.649596968324358   1.395449806257184   2.332526365071812   2.150211932758257];
P.phis2 = x(1); % increasing function in immunity
P.phir2 = x(2); 
P.rhos2 = x(3);
P.rhor2 = x(4); 
P.psis2 = x(5);
P.psir2 = x(6);

% % dynamical immunity
P.phif0 = 0.01; P.phif1 = 1;
P.rhof0 = 0.01; P.rhof1 = 1; 
P.psif0 = 0.01; P.psif1 = 1; 

% % Fixed immunity (max=min=constant)
% P.phif0 = 0.5; P.phif1 = 0.5; 
% P.rhof0 = 0.5; P.rhof1 = 0.5; 
% P.psif0 = 0.5; P.psif1 = 0.5; 

% Fixed high immunity (max=min=constant)
% P.phif0 = 0.9; P.phif1 = 0.9; 
% P.rhof0 = 0.9; P.rhof1 = 0.9; 
% P.psif0 = 0.9; P.psif1 = 0.9; 

% Fixed low immunity (max=min=constant)
% P.phif0 = 0.2; P.phif1 = 0.2; 
% P.rhof0 = 0.2; P.rhof1 = 0.2; 
% P.psif0 = 0.2; P.psif1 = 0.2; 

[P.rho, P.phi, P.psi] = sigmoid_prob(0, P);

%%
% figure_setups; hold on
% xx = 1:0.01:15;
% [rho, phi, psi] = sigmoid_prob(xx, P);
% plot(xx,rho)
% plot(xx,phi)
% plot(xx,psi)
% legend('$\rho$', '$\phi$','$\psi$')
