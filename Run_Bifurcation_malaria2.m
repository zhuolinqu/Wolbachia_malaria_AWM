%% bifurcation plot (R0M vs. malaria prevalence)
clearvars; 
% close all; clc
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
R0M = NaN(6,length(phiW_list));
for iphi = 1:length(phiW_list)   
    P.phiW = phiW_list(iphi);
    if iphi == 1
        SS_mat = EquilibriumState_m(P);
    else
        SS_mat = EquilibriumState_m(P,SS_mat_old);
    end
    Minf(:,iphi) = (SS_mat(:,3)+SS_mat(:,4))./sum(SS_mat(:,1:4),2);
    Stab(:,iphi,:) = SS_mat(:,12:13);
    SS_matW = EquilibriumState_w(P);
    R0M(1,iphi) = Cal_R0_malaria(SS_matW(1,1),SS_matW(1,2),P);
    R0M(2,iphi) = Cal_R0_malaria(SS_matW(2,1),SS_matW(2,2),P);
    R0M(3,iphi) = Cal_R0_malaria(SS_matW(3,1),SS_matW(3,2),P);
    R0M(4,iphi) = R0M(1,iphi);
    R0M(5,iphi) = R0M(2,iphi);
    R0M(6,iphi) = R0M(3,iphi);
    SS_mat_old = SS_mat;
end
toc
legend_list = {'stable','Wolbachia-unstable (malaria stable)','malaria-unstable (Wolbachia-stable)',...
    'malaria-unstable \& Wolbachia-unstable'};
%%
f = figure_setups;
hold on
for iline = 1:6
    group1 = find((Stab(iline,:,1)==1).*(Stab(iline,:,2)==1));
    group2 = find((Stab(iline,:,1)==1).*(1-Stab(iline,:,2)==1));
    group3 = find((1-Stab(iline,:,1)==1).*(Stab(iline,:,2)==1));
    group4 = find((1-Stab(iline,:,1)==1).*(1-Stab(iline,:,2)==1));   
    plot(R0M(iline, group1),Minf(iline,group1),'-','Color',[0 0.4470 0.7410],'DisplayName',legend_list{1})
    plot(R0M(iline, group2),Minf(iline,group2),'-.','Color',[0.4660 0.6740 0.1880],'DisplayName',legend_list{2})
    plot(R0M(iline, group3),Minf(iline,group3),'--','Color',[0.8500 0.3250 0.0980],'DisplayName',legend_list{3})
    plot(R0M(iline, group4),Minf(iline,group4),':','Color',[0.6350 0.0780 0.1840],'DisplayName',legend_list{4})
end
ll = legendUnq(f);
ll = ll([3,4,1,2]);
legend(ll,'Location','nw')
xlabel('$\mathcal{R}_0^m$')
ylabel('Malaria prevalence')
title(['$v_w=',num2str(P.vw),',~~ c_i=', num2str(P.ci),'$'])
xlim([0.5 3])

