%% bifurcation plot (R0M vs. varying R0W)
clearvars;
close all; clc
tic

%% Parameters & numerical config
Baseline_params_malaria;
P = Baseline_params_stephensi(P);

for izone = 1:6
    label_zone = izone;
    [P.vw,P.ci] = param_select(label_zone); P.vu = 1- P.vw;
    % R0w list
    R0w_del = Cal_R0w_Del(P);
    R0w_max = 2;
    R0W_list = [linspace(0.01,R0w_del-10^-3,50),linspace(R0w_del+10^-3,0.99,50),linspace(1.01,R0w_max,50)];
    %% Run steady state calculations
    R0M_list = NaN(3,length(R0W_list));
    Stab_list = NaN(3,length(R0W_list));
    for iR0W = 1:length(R0W_list)
        R0W = R0W_list(iR0W);
        % P.phiU = P.vw*P.phiW*P.mufu/(P.mufw*R0W);
        P.phiW = (P.mufw*R0W*P.phiU)/(P.vw*P.mufu);
        SS_mat = EquilibriumState_w(P);
        Stab_list(:,iR0W) = SS_mat(:,end);
        [R0M_list(1,iR0W), ~, ~] = Cal_R0_malaria(SS_mat(1,1),SS_mat(1,2),P);
        [R0M_list(2,iR0W), ~, ~] = Cal_R0_malaria(SS_mat(2,1),SS_mat(2,2),P);
        [R0M_list(3,iR0W), ~, ~] = Cal_R0_malaria(SS_mat(3,1),SS_mat(3,2),P);
    end
    toc
    %
    figure_setups; hold on
    for i=1:3
        ind_stab = find(Stab_list(i,:)==1);
        ind_unstab = find(Stab_list(i,:)==0);
        plot(R0W_list(ind_stab),R0M_list(i,ind_stab),'k-')
        plot(R0W_list(ind_unstab),R0M_list(i,ind_unstab),'k--')
    end
    plot([0 2],[1 1],'r:')
    plot([1 1],[0 5],'r:')
    xlabel('$\mathcal{R}_0^w$')
    ylabel('$\mathcal{R}_0^m$')
    axis([0 R0w_max 0 5])
    title(['Zone ',num2str(label_zone)])
end
%%
function [vw,ci] = param_select(label_zone)
switch label_zone
    case 1
        ci = 0.5; vw = 0.6;
    case 2
        ci = 0.9; vw = 0.9;
    case 3
        ci = 0.01; vw = 0.98;
    case 4
        ci = 0.2; vw = 0.98;
    case 5
        ci = 0.95; vw = 0.96;
    case 6
        ci = 0.97; vw = 0.99;
end
end
