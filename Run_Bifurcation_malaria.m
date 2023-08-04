clearvars; 
% close all; clc
tic

%% Parameters & numerical config
Baseline_params_malaria;
P = Baseline_params_stephensi(P);
P.vw = 0.95; P.vu = 1- P.vw;
%% Sampling
phiW_list = linspace(0.8,0.85,1000);

%% Run steady state calculations
Minf = NaN(6,length(phiW_list));
Minf_stab = NaN(6,length(phiW_list));
R0w = NaN(1,length(phiW_list));
for iphi = 1:length(phiW_list)   
    P.phiW = phiW_list(iphi);
    SS_mat = EquilibriumState_m(P);
    Minf(:,iphi) = (SS_mat(:,3)+SS_mat(:,4))./sum(SS_mat(:,1:4),2);
    Minf_stab(:,iphi) = SS_mat(:,end);
    R0w(:,iphi) = Cal_R0_wolbachia(P);
end

legend_list = {'no malaria and no Wol.','no malaria and unstable Wol.','no malaria and stable Wol',...
    'malaria endemic and no Wol','malaria endemic and unstable Wol','malaria endemic and stable Wol'};

figure
hold on
for iline = 1:6
    legend_list{iline}
    disp('stable')
    plot(R0w(Minf_stab(iline,:)==1),Minf(iline,Minf_stab(iline,:)==1),'x-','DisplayName',legend_list{iline},'LineWidth',2)
    pause
    disp('unstable')
    plot(R0w(Minf_stab(iline,:)==0),Minf(iline,Minf_stab(iline,:)==0),'--','DisplayName',legend_list{iline},'LineWidth',2)
    pause
end
legend
xlabel('R0w')
ylabel('Malaria prevalence')

toc
