%% bifurcation plot (R0w vs. wolbachia prevalence)
clearvars; 
% close all; clc
tic

%% Parameters & numerical config
Baseline_params_malaria;
P = Baseline_params_stephensi(P);
P.vw = 0.95; P.vu = 1- P.vw;
P.ci = 1;
%% Sampling
phiW_min = P.mufw/(P.vw*P.bf);
phiW_list = [phiW_min:0.001:0.245, 0.245:0.001:0.25, 0.25:0.005:0.5, ...
    0.5:0.005:0.55, 0.55:0.01:1.5]; 
% phiW_list = linspace(0.01,2,10000);

%% Run model
Winf = NaN(3,length(phiW_list));
Winf_stab = NaN(3,length(phiW_list));
R0w = NaN(1,length(phiW_list));

for iphi = 1:length(phiW_list)   
    P.phiW = phiW_list(iphi);
    SS_mat = EquilibriumState_w(P);
    Winf(:,iphi) = SS_mat(:,2)./sum(SS_mat(:,1:end-1),2);
    Winf_stab(:,iphi) = SS_mat(:,end);
    R0w(:,iphi) = Cal_R0_wolbachia(P);
end

if P.vw<1
    legend_list = {'WFE','WEE$^-$','WEE$^+$'};
elseif P.vw==1
    legend_list = {'WFE','WEE$^-$','WCE'};
end

figure_setups
hold on

for iline = 1:3  
    plot(R0w(Winf_stab(iline,:)==1),Winf(iline,Winf_stab(iline,:)==1),'-','DisplayName',legend_list{iline})   
    plot(R0w(Winf_stab(iline,:)==0),Winf(iline,Winf_stab(iline,:)==0),'--','DisplayName',legend_list{iline})
end
legend
xlabel('$\mathcal{R}_0^w$')
ylabel('Wolbachia prevalence')
axis([0 1.2 0 1.1])
title(['$v_w=',num2str(P.vw),',~~ c_i=', num2str(P.ci),'$'])
toc