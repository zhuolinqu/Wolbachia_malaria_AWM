clearvars; 
% close all; clc
tic

%% Parameters & numerical config
Baseline_params_malaria;
P = Baseline_params_stephensi(P);
P.vw = 0.95; P.vu = 1- P.vw;
%% Sampling
phiW_list = linspace(0.1,2,50);%  [linspace(0.2,0.85,100), linspace(0.8,0.85,200),linspace(0.85,5,200)];

%% Run steady state calculations
Minf = NaN(6,length(phiW_list));
Minf_stab = NaN(6,length(phiW_list));
R0w = NaN(1,length(phiW_list));
for iphi = 1:length(phiW_list)   
    P.phiW = phiW_list(iphi);
    if iphi == 1
        SS_mat = EquilibriumState_m(P);
    else
        SS_mat = EquilibriumState_m(P,SS_mat_old);
    end
    Minf(:,iphi) = (SS_mat(:,3)+SS_mat(:,4))./sum(SS_mat(:,1:4),2);
    Minf_stab(:,iphi) = SS_mat(:,end);
    R0w(:,iphi) = Cal_R0_wolbachia(P);
    SS_mat_old = SS_mat;
end

legend_list = {'no malaria and no Wol.','no malaria and unstable Wol.','no malaria and stable Wol',...
    'malaria endemic and no Wol','malaria endemic and unstable Wol','malaria endemic and stable Wol'};

figure_setups

hold on
for iline = 1:6
    legend_list{iline}
    disp('stable')
    plot(R0w(Minf_stab(iline,:)==1),Minf(iline,Minf_stab(iline,:)==1),'x-','DisplayName',legend_list{iline})
    % pause
    disp('unstable')
    plot(R0w(Minf_stab(iline,:)==0),Minf(iline,Minf_stab(iline,:)==0),'--','DisplayName',legend_list{iline})
    % pause
end
% legend
xlabel('$R_0^w$')
ylabel('Malaria prevalence')

toc
