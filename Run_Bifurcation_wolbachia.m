%% bifurcation plot (R0w vs. wolbachia prevalence)
clearvars; 
close all; clc
tic

%% Parameters & numerical config
Baseline_params_malaria;
P = Baseline_params_stephensi(P);
P.vw = 0.95; P.vu = 1- P.vw;
P.ci = 0.9;
% P.vw = 1; P.vu = 1- P.vw;
% P.ci = 0.8;
%% Sampling
phiW_min = P.mufw/(P.vw*P.bf);
phiW_list = [0.01:0.001:0.245, 0.245:0.001:0.25, 0.25:0.005:0.5, ...
    0.5:0.005:0.55, 0.55:0.01:4]; 

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

legend_list = {'W-SS stable','W-SS unstable'};
%%
f = figure_setups;
hold on
for iline = 1:3  
    group1 = find((Winf_stab(iline,:,1)==1));
    group2 = find((Winf_stab(iline,:,1)==0));
    plot(R0w(Winf_stab(iline,:)==1),Winf(iline,Winf_stab(iline,:)==1),'-','Color',[0 0.4470 0.7410],'DisplayName',legend_list{1})   
    plot(R0w(Winf_stab(iline,:)==0),Winf(iline,Winf_stab(iline,:)==0),'-.','Color',[0.4660 0.6740 0.1880],'DisplayName',legend_list{2})
end
ll = legendUnq(f);
legend(ll,'Position',[0.15,0.34,0.32,0.13])
xlabel('$\mathcal{R}_0^w$')
ylabel('Wolbachia prevalence')
axis([0 1.2 0 1.1])

if P.vw<1
    text(0.25,0.05,'WFE')
    text(0.8,0.3,'WEE$^-$') 
    text(0.8,0.87,'WEE$^+$') 
elseif P.vw==1
    text(0.25,0.05,'WFE')
    text(0.8,0.4,'WEE$^-$')  
    text(0.8,0.9,'WCE') 
    % legend_list = {'WFE','WEE$^-$','WCE'};
end

title(['$v_w=',num2str(P.vw),',~~ c_i=', num2str(P.ci),'$'])
toc