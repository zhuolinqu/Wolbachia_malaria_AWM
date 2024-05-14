%% bifurcation plot (R0M vs. malaria prevalence) varying betaM
clearvars; 
close all; clc
tic

%% Parameters & numerical config
Baseline_params_malaria;
P = Baseline_params_stephensi(P);
P.vw = 0.95; P.vu = 1- P.vw;
Cal_R0_wolbachia(P)
SS_matW = EquilibriumState_w(P); SU = SS_matW(1,1);
Cal_R0_malaria(SU,0,P)

[betaM_min, R0m_min] = Cal_betaM_R0m(P,0.5);
[betaM1, ~] = Cal_betaM_R0m(P,1);
[betaM2, ~] = Cal_betaM_R0m(P,1.5);
[betaM_max, R0m_max] = Cal_betaM_R0m(P,3);
%% Sampling
betaM_list = [linspace(betaM_min,betaM1,3),linspace(betaM1,betaM2,10),linspace(betaM2,betaM_max,20)];

%% Run steady state calculations
Minf = NaN(1,length(betaM_list));
R0M = NaN(1,length(betaM_list));
for ibetaM = 1:length(betaM_list)   
    P.betaM = betaM_list(ibetaM);
    [R0m, SS_mat] = EquilibriumState_m_row4(P);
    Minf(1,ibetaM) = (SS_mat(1,3)+SS_mat(1,4))./sum(SS_mat(1,1:4),2);
    R0M(1,ibetaM) = R0m;
end
toc
%%
figure_setups;
hold on
ind = find(R0M>=0.99);
plot(R0M(ind),Minf(ind),'-','DisplayName','malaria stable')
ind = find(R0M<1);
plot(R0M(ind),Minf(ind),'--','DisplayName','malaria unstable')
xlabel('$\mathcal{R}_0^m$')
ylabel('Malaria prevalence')
xlim([R0m_min R0m_max])
legend('Location','nw')

function [R0m, SS_mat] = EquilibriumState_m_row4(P)
[R0w, G0w, G0u] = Cal_R0_wolbachia(P);

if G0u<1 || G0w<1
    disp('mosquito extinction')
    % mosquito extinction, check parameter range!
end

SS_mat = NaN(1,11);

SS_matW = EquilibriumState_w(P);
Wol_DFE = SS_matW(1,1:end-1); 
SU = Wol_DFE(1); SW = Wol_DFE(2); % SW = 0;
R0m = Cal_R0_malaria(SU,SW,P);

SH0 = P.gH/P.muH-1; EH0 = 1; AH0 = 0; DH0 = 0; Ie0 = 0;
SU0 = SU; EU0 = 0; IU0 = 0; SW0 = 0; EW0 = 0; IW0 = 0;
yinit = [SH0; EH0; AH0; DH0; Ie0; SU0; EU0; IU0; SW0; EW0; IW0];
options = odeset('AbsTol',1e-10,'RelTol',1e-10);
[~,y] = ode45(@BaseModel,linspace(0,10^4,50),yinit,options,P);
SS_mat(1,1:11) = y(end,:);

end