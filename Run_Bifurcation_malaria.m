clearvars; 
% close all; clc
tic

%% Parameters & numerical config
Baseline_params_malaria;
P = Baseline_params_stephensi(P);
P.vw = 0.95; P.vu = 1- P.vw;
%% Sampling
phiW_list = linspace(0.2,5,500);

%% Run model
Minf = NaN(6,length(phiW_list));
R0w = NaN(1,length(phiW_list));
for iphi = 1:length(phiW_list)   
    P.phiW = phiW_list(iphi);
    SS_mat = EquilibriumState_m(P);
    Minf(:,iphi) = (SS_mat(:,3)+SS_mat(:,4))./sum(SS_mat(:,1:4),2);
    R0w(:,iphi) = Cal_R0_wolbachia(P);
end

figure
plot(R0w,Minf,'*','linewidth',2)
xlabel('R0w')
ylabel('Malaria prevalence')
legend('no malaria and no Wol.','no malaria and unstable Wol.','no malaria and stable Wol',...
    'malaria endemic and no Wol','malaria endemic and unstable Wol','malaria endemic and stable Wol')
toc

% (row 1) DFE-DFE: no malaria and no Wolbachia
% (row 2) DFE-EE-: no malaria and unstable Wolbachia endemic
% (row 3) DFE-EE+: no malaria and stable Wolbachia endemic
% (row 4) EE-DFE: malaria endemic and no Wolbachia
% (row 5) EE-EE-: malaria endemic and unstable Wolbachia endemic
% (row 6) EE-EE+: malaria endemic and stable Wolbachia endemic
