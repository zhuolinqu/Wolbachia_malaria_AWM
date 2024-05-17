%% bifurcation plot (R0w vs. malaria prevalence)
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

%% Sampling
phiW_min = P.mufw/(P.vw*P.bf);
phiW_list = [phiW_min:0.03:0.415, 0.415:0.001:0.425, ...
    0.425:0.005:0.55, 0.55:0.05:3];  

%% Run steady state calculations
Minf = NaN(6,length(phiW_list));
Stab = NaN(6,length(phiW_list),2);
R0w = NaN(1,length(phiW_list));
for iphi = 1:length(phiW_list)   
    P.phiW = phiW_list(iphi);
    if iphi == 1
        SS_mat = EquilibriumState_m(P);
    else
        SS_mat = EquilibriumState_m(P,SS_mat_old);
    end
    Minf(:,iphi) = (SS_mat(:,3)+SS_mat(:,4))./sum(SS_mat(:,1:4),2);
    Stab(:,iphi,:) = SS_mat(:,12:13);
    R0w(:,iphi) = Cal_R0_wolbachia(P);
    SS_mat_old = SS_mat;
end
toc
legend_list = {'W-SS stable, malaria-SS stable','W-SS unstable, malaria-SS stable','W-SS stable, malaria-SS unstable',...
    'W-SS unstable, malaria-SS unstable'};
%%
f = figure_setups;
hold on
for iline = 1:6
    group1 = find((Stab(iline,:,1)==1).*(Stab(iline,:,2)==1));
    group2 = find((Stab(iline,:,1)==1).*(1-Stab(iline,:,2)==1));
    group3 = find((1-Stab(iline,:,1)==1).*(Stab(iline,:,2)==1));
    group4 = find((1-Stab(iline,:,1)==1).*(1-Stab(iline,:,2)==1));   
    plot(R0w(group1),Minf(iline,group1),'-','Color',[0 0.4470 0.7410],'DisplayName',legend_list{1})
    plot(R0w(group2),Minf(iline,group2),'-.','Color',[0.4660 0.6740 0.1880],'DisplayName',legend_list{2})
    plot(R0w(group3),Minf(iline,group3),'--','Color',[0.8500 0.3250 0.0980],'DisplayName',legend_list{3})
    plot(R0w(group4),Minf(iline,group4),':','Color',[0.6350 0.0780 0.1840],'DisplayName',legend_list{4})
end
ll = legendUnq(f);
ll = ll([3,4,1,2]);
legend(ll,'Location','east')
xlabel('$\mathcal{R}_0^w$')
ylabel('Malaria prevalence')
ylim([0 0.8])
print(gcf,'-vector', '-depsc', 'Results/bifur_fixed.eps')

% print(gcf,'-vector', '-depsc', 'Results/bifur_dynamic.eps')
