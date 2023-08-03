clearvars; close all; clc
tic

%% Parameters & numerical config
Baseline_params_malaria;
P = Baseline_params_stephensi(P);
P.vw = 1; P.vu = 1- P.vw;
%% Sampling
phiW_list = linspace(0,10,100000);

%% Run model
Winf = NaN(3,length(phiW_list));
R0w = NaN(1,length(phiW_list));
for iphi = 1:length(phiW_list)   
    P.phiW = phiW_list(iphi);
    SS_mat = EquilibriumState_w(P);
    Winf(:,iphi) = SS_mat(:,2)./sum(SS_mat,2);
    R0w(:,iphi) = Cal_R0_wolbachia(P);
    
end

figure(1)
plot(R0w,Winf,'-','linewidth',2)
toc

