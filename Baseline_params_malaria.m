%% Malaria parameters taken from Qu and Patterson et al 2023 (except as noted)
P.bm = 0.7;
P.bh = 5;

P.betaM = 0.3; P.betaM_lower = 0.21; P.betaM_upper = 0.51;
P.betaD = 0.35; P.betaD_lower = 0.25; P.betaD_upper = 0.6;
P.betaA = 0.03; P.betaA_lower = 0.02; P.betaA_upper = 0.05;

P.rD = 1/33.5; % Updated based on treatment assumption using malaria therapy data % 1/180; original Qu and Patterson et al
P.rD_lower = 1/48; P.rD_upper = 1/20;
P.rA = 1/85; % Updated  using malaria therapy data  % 1/360; original Qu and Patterson et al
P.rA_lower = 1/121; P.rA_upper = 1/50;

P.h = 1/26; % Updated based on malaria therapy data % 1/15; original Qu and Patterson et al
P.h_lower = 1/37; P.h_upper = 1/15;
P.muD = 0;

P.de = 5*365; P.de_lower = 3.25*365; P.de_upper = 10.5*365;
P.gamma = 10;

% Simple birth/death parameters (ave lifespan = 60)
P.gH = 1;
P.muH = 1/(65*365);

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
% P.phif0 = 0.05; P.phif1 = 0.05; 
% P.rhof0 = 0.05; P.rhof1 = 0.05; 
% P.psif0 = 0.05; P.psif1 = 0.05; 

% sigmod parameters for comparison
% P.phif0 = 0.329292647824602; P.phif1 = 0.329292647824602; 
% P.rhof0 = 0.118177744193668; P.rhof1 = 0.118177744193668; 
% P.psif0 = 0.327791291717994; P.psif1 = 0.327791291717994; 

[P.rho, P.phi, P.psi] = sigmoid_prob(0, P);

low = 0.7; high = 1.7;

P.phis2_lower = P.phis2*low; P.phis2_upper = P.phis2*high;
P.phir2_lower = P.phir2*low; P.phir2_upper = P.phir2*high;
P.rhos2_lower = P.rhos2*low; P.rhos2_upper = P.rhos2*high;
P.rhor2_lower = P.rhor2*low; P.rhor2_upper = P.rhor2*high;
P.psis2_lower = P.psis2*low; P.psis2_upper = P.psis2*high;
P.psir2_lower = P.psir2*low; P.psir2_upper = P.psir2*high;

%%
P.cS = 0.75;
P.cE = 0.1;
P.cA = 0.1;
P.cD = 0.05;

P.cS_lower = P.cS*low; P.cS_upper = P.cS*high;
P.cE_lower = P.cE*low; P.cE_upper = P.cE*high;
P.cA_lower = P.cA*low; P.cA_upper = P.cA*high;
P.cD_lower = P.cD*low; P.cD_upper = P.cD*high;

%% dummy varaible for eFast SA
P.dummy = 1; P.dummy_lower = 0.7; P.dummy_upper = 1.7; % dummy parameter for global SA eFAST - values from original code 

%%
% figure_setups; hold on
% xx = 1:0.01:15;
% [rho, phi, psi] = sigmoid_prob(xx, P);
% plot(xx,rho)
% plot(xx,phi)
% plot(xx,psi)
% legend('$\rho$', '$\phi$','$\psi$')
