%% bifurcation plot (R0w vs. malaria prevalence vs. wolbachia prevalence)
clearvars;
% close all; clc
tic
%% Parameters & numerical config
Baseline_params_malaria;
P = Baseline_params_stephensi(P);
P.vw = 0.95; P.vu = 1- P.vw;
%% Sampling
phiW_min = P.mufw/(P.vw*P.bf);
phiW_list = [phiW_min:0.03:0.145, 0.145:0.0001:0.150, 0.15:0.005:0.5, ...
    0.5:0.005:0.55, 0.55:0.01:1.5]; 

%% Run steady state calculations
tic
Minf = NaN(6,length(phiW_list));
Winf = NaN(6,length(phiW_list));
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
    Winf(:,iphi) = sum(SS_mat(:,9:11),2)./sum(SS_mat(:,6:11),2);
    Stab(:,iphi,:) = SS_mat(:,12:13);
    R0w(:,iphi) = Cal_R0_wolbachia(P);
    SS_mat_old = SS_mat;
end
save('Bifurcation_3D_low.mat')
toc
%% plotting
% load('Bifurcation_3D.mat')
f = figure_setups; hold on;
set(f,'WindowState','maximized') % make the plot full screen
view([65,45])
axis([0 1 0 1.2 0 1])
xlabel('Wolbachia prevalence')
ylabel('$R_0^w$')
zlabel('Malaria prevalence $(A_H+D_H)$')
grid on
legend_list = {'stable','Wolbachia-unstable (malaria stable)','malaria-unstable (Wolbachia-stable)',...
    'malaria-unstable \& Wolbachia-unstable'};
for iline = 1:6
    group1 = find((Stab(iline,:,1)==1).*(Stab(iline,:,2)==1));
    group2 = find((Stab(iline,:,1)==1).*(1-Stab(iline,:,2)==1));
    group3 = find((1-Stab(iline,:,1)==1).*(Stab(iline,:,2)==1));
    group4 = find((1-Stab(iline,:,1)==1).*(1-Stab(iline,:,2)==1));   
    plot3(Winf(iline,group1),R0w(group1),Minf(iline,group1),'-','Color',[0 0.4470 0.7410],'DisplayName',legend_list{1})
    plot3(Winf(iline,group2),R0w(group2),Minf(iline,group2),'-.','Color',[0.4660 0.6740 0.1880],'DisplayName',legend_list{2})
    plot3(Winf(iline,group3),R0w(group3),Minf(iline,group3),'--','Color',[0.8500 0.3250 0.0980],'DisplayName',legend_list{3})
    plot3(Winf(iline,group4),R0w(group4),Minf(iline,group4),':','Color',[0.6350 0.0780 0.1840],'DisplayName',legend_list{4})
end
ll = legendUnq(f);
ll = ll([3,4,1,2]);
legend(ll,'Location','east')

%% 
% legend_list = {'no malaria and no Wol.','no malaria and unstable Wol.','no malaria and stable Wol',...
%     'malaria endemic and no Wol','malaria endemic and unstable Wol','malaria endemic and stable Wol'};
% for iline = 1:6
%     plot3(Winf(iline,Stab(iline,:)==1),R0w(Stab(iline,:)==1),Minf(iline,Stab(iline,:)==1),'-','DisplayName',legend_list{iline})
%     plot3(Winf(iline,Stab(iline,:)==0),R0w(Stab(iline,:)==0),Minf(iline,Stab(iline,:)==0),'--','DisplayName',legend_list{iline})
% end
% legend

% %% solution trajectory
% P = Baseline_params_stephensi(P);
% phi_list = 0.1:0.2:0.9;
% for phi = phi_list    
%     P.vw = 0.95; P.vu = 1- P.vw;
%     P.phiW = phi;
%     SS_mat = EquilibriumState_m(P);
%     R0w = Cal_R0_wolbachia(P);
%     yinit = SS_mat(5,1:11); 
%     %[SH, EH, AH, DH, Ie, SU, EU, IU, SW, EW, IW]
%     if isnan(yinit)
%         yinit = SS_mat(1,1:11);
%         yinit(2) = 0.1*yinit(1); % EH
%         yinit(1) = 0; % SH
%         yinit(9) = yinit(6)*0.9; % SW
%         yinit(6) = yinit(6)*0.1; % SU
%     else
%         yinit(1) = yinit(1)+yinit(3)+yinit(4); % SH
%         yinit(3) = 0; % AH
%         yinit(4) = 0; % DH
%         yinit(9) = yinit(9)-yinit(6)*0.8; % SW
%         yinit(6) = yinit(6)+yinit(6)*0.8; % SU
%     end
%     options = odeset('AbsTol',1e-10,'RelTol',1e-10);
%     [t,y] = ode45(@BaseModel,linspace(0,3000,1500),yinit,options,P);
% 
%     SH = y(:,1); EH = y(:,2); AH = y(:,3); DH = y(:,4); Ie = y(:,5); 
%     SU = y(:,6); EU = y(:,7); IU = y(:,8); SW = y(:,9); EW = y(:,10); IW = y(:,11);
%     NH = SH + EH + AH + DH; NU = SU + EU + IU; NW = SW + EW + IW; NM = NU + NW;
% 
%     col = t;
%     scatter3(NW./NM,R0w*ones(size(t)),(AH+DH)./NH,[],col,'>','filled')
% end
% %%
% phi_list = [1.3, 1, 0.7];
% for phi = phi_list    
%     P.vw = 0.95; P.vu = 1- P.vw;
%     P.phiW = phi;
%     SS_mat = EquilibriumState_m(P);
%     R0w = Cal_R0_wolbachia(P);
%     yinit = SS_mat(5,1:11); 
%     %[SH, EH, AH, DH, Ie, SU, EU, IU, SW, EW, IW]
%     if isnan(yinit)
%         yinit = SS_mat(1,1:11);
%         yinit(2) = 0.1*yinit(1); % EH
%         yinit(1) = 0; % SH
%         yinit(9) = yinit(6)*0.1; % SW
%         yinit(6) = yinit(6)*0.9; % SU
%     else
%         yinit(1) = yinit(1)+yinit(3)+yinit(4); % SH
%         yinit(3) = 0; % AH
%         yinit(4) = 0; % DH
%         yinit(9) = yinit(9)+yinit(6)*0.2; % SW
%         yinit(6) = yinit(6)-yinit(6)*0.2; % SU
%     end
%     [t,y] = ode45(@BaseModel,linspace(0,3000,1500),yinit,options,P);
% 
%     SH = y(:,1); EH = y(:,2); AH = y(:,3); DH = y(:,4); Ie = y(:,5); 
%     SU = y(:,6); EU = y(:,7); IU = y(:,8); SW = y(:,9); EW = y(:,10); IW = y(:,11);
%     NH = SH + EH + AH + DH; NU = SU + EU + IU; NW = SW + EW + IW; NM = NU + NW;
% 
%     col = t;
%     scatter3(NW./NM,R0w*ones(size(t)),(AH+DH)./NH,[],col,'>','filled')
% end
