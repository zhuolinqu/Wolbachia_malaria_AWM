clearvars;
close all; clc
tic


%% Parameters & numerical config
% Baseline_params_malaria;
% P = Baseline_params_stephensi(P);
% P.vw = 0.95; P.vu = 1- P.vw;
% %% Sampling
% phiW_min = P.mufw/(P.vw*P.bf);
% phiW_list = [phiW_min:0.03:0.245, 0.245:0.001:0.25, 0.25:0.005:0.5, ...
%     0.5:0.005:0.55, 0.55:0.01:1.5];   %  [linspace(0.2,0.85,100), linspace(0.8,0.85,200),linspace(0.85,5,200)]; 
% 
% %% Run steady state calculations
% Minf = NaN(6,length(phiW_list));
% Winf = NaN(6,length(phiW_list));
% Stab = NaN(6,length(phiW_list));
% R0w = NaN(1,length(phiW_list));
% for iphi = 1:length(phiW_list)
%     P.phiW = phiW_list(iphi);
%     SS_mat = EquilibriumState_m(P);
%     Minf(:,iphi) = (SS_mat(:,3)+SS_mat(:,4))./sum(SS_mat(:,1:4),2);
%     Winf(:,iphi) = sum(SS_mat(:,9:11),2)./sum(SS_mat(:,6:11),2);
%     Stab(:,iphi) = SS_mat(:,end);
%     R0w(:,iphi) = Cal_R0_wolbachia(P);
% end
% save('Bifurcation_3D.mat')

load('Bifurcation_3D.mat')
%%
% legend_list = {'no malaria and no Wol.','no malaria and unstable Wol.','no malaria and stable Wol',...
%     'malaria endemic and no Wol','malaria endemic and unstable Wol','malaria endemic and stable Wol'};
% 
% figure_setups
% hold on
% for iline = 1:6
%     plot3(R0w(Stab(iline,:)==1),Minf(iline,Stab(iline,:)==1),Winf(iline,Stab(iline,:)==1),'x-','DisplayName',legend_list{iline},'LineWidth',2)
%     plot3(R0w(Stab(iline,:)==0),Minf(iline,Stab(iline,:)==0),Winf(iline,Stab(iline,:)==0),'--','DisplayName',legend_list{iline},'LineWidth',2)
% end
% legend
% xlabel('R0w')
% ylabel('Malaria prevalence')
% zlabel('Wolbachia prevalence')
% view(3)
% grid on
% toc

%% change angle to match scratch paper
legend_list = {'no malaria and no Wol.','no malaria and unstable Wol.','no malaria and stable Wol',...
    'malaria endemic and no Wol','malaria endemic and unstable Wol','malaria endemic and stable Wol'};

f = figure_setups; hold on; 
set(f,'WindowState','maximized') % make the plot full screen
view([110,17]) 

for iline = 1:6
    plot3(Winf(iline,Stab(iline,:)==1),R0w(Stab(iline,:)==1),Minf(iline,Stab(iline,:)==1),'-','DisplayName',legend_list{iline})
    plot3(Winf(iline,Stab(iline,:)==0),R0w(Stab(iline,:)==0),Minf(iline,Stab(iline,:)==0),'--','DisplayName',legend_list{iline})
end
% legend
xlabel('Wolbachia prevalence')
ylabel('$R_0^w$')
zlabel('Malaria prevalence')
grid on

%% solution trajectory
% rest to baseline parameters
P = Baseline_params_stephensi(P); 
P.vw = 0.95; P.vu = 1- P.vw;
SS_mat = EquilibriumState_m(P);
R0w = Cal_R0_wolbachia(P);
yinit = SS_mat(5,1:end-1);
%[SH, EH, AH, DH, Ie, SU, EU, IU, SW, EW, IW]
yinit(1) = yinit(1)+yinit(3)+yinit(4)+yinit(2)*0;
yinit(2) = yinit(2)-yinit(2)*0;
yinit(3) = yinit(3)-yinit(3);
yinit(4) = yinit(4)-yinit(4);
yinit(9) = yinit(9)+yinit(6)*0.8; % perturb the wolbachia
yinit(6) = yinit(6)-yinit(6)*0.8;
options = odeset('AbsTol',1e-10,'RelTol',1e-10);
[t,y] = ode45(@BaseModel,linspace(0,3000,1500),yinit,options,P);

SH = y(:,1);
EH = y(:,2);
AH = y(:,3);
DH = y(:,4);
Ie = y(:,5);
SU = y(:,6);
EU = y(:,7);
IU = y(:,8);
SW = y(:,9);
EW = y(:,10);
IW = y(:,11);

NH = SH + EH + AH + DH;
NU = SU + EU + IU;
NW = SW + EW + IW;
NM = NU + NW;

col = t;
scatter3(NW./NM,R0w*ones(size(t)),(AH+DH)./NH,[],col,'>','filled')
% colorbar
% caxis([0 5000])


