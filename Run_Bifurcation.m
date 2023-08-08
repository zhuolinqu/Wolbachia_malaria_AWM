clearvars; close all; clc
%global P

tic

%% Parameters & numerical config
Baseline_params_malaria;
P = Baseline_params_stephensi(P);
tfinal = 1000;
%% Sampling
EH0_max = P.gH/P.muH; EH0_min = 1;
SW0_max = 0.8*P.Kf; SW0_min = 1;

% LHS_raw = lhsdesign(NS,3);
%
% EH0_list = LHS_raw(:,1)*(EH0_max-EH0_min)+EH0_min;
% SW0_list = LHS_raw(:,2)*(SW0_max-SW0_min)+SW0_min;
% phiW_list = LHS_raw(:,3)*(phiU_max-phiU_min)+phiU_min;
phiW_list = linspace(5,15,10);
EH0_list = linspace(EH0_min,EH0_max,10);
SW0_list = linspace(SW0_min,SW0_max,10);

%% Run model
% Initial conditions
AH0 = 0; DH0 = 0; Ie0 = 0;
EU0 = 0; IU0 = 0; EW0 = 0; IW0 = 0;

for iphi = 1:length(phiW_list)   
    for iEH = 1:length(EH0_list)
        for iSW = 1:length(SW0_list)
            [iphi,iEH,iSW]
            
            P.phiW = phiW_list(iphi);
            EH0 = EH0_list(iEH); SH0 = P.gH/P.muH-EH0;
            SW0 = SW0_list(iSW); SU0 = P.Kf-SW0;
            
            yinit = [SH0; EH0; AH0; DH0; Ie0;
                SU0; EU0; IU0; ...
                SW0; EW0; IW0];
            
            R0_wol = Cal_R0_wolbachia(P);
            
            [t,y] = ode45(@BaseModel,linspace(0,tfinal,100),yinit,[],P);
            t_cut = 15:length(t);
            SH = y(:,1); EH = y(:,2); AH = y(:,3); DH = y(:,4);
            SU = y(:,6); EU = y(:,7); IU = y(:,8);
            SW = y(:,9); EW = y(:,10); IW = y(:,11);
            p_w = (SW+EW+IW)./sum(y(:,6:11),2); % wolbachia prevalence
            p_m = (AH+DH)./sum(y(:,1:4),2); % malaria prevalence
            % remove the initial transient dynamics
            t = t(t_cut);
            p_w = p_w(t_cut);
            p_m = p_m(t_cut);
            
            % plot
            if p_w(end)>10^-2; lcolor = 'r'; else; lcolor = 'b'; end
            figure(1);
            subplot(1,2,1); hold on
            scatter(R0_wol*ones(size(p_w)),p_w,20,'MarkerEdgeColor',lcolor,'MarkerFaceColor',lcolor)
            hold off
            ylim([0 1])
            xlabel('R0 wol')
            ylabel('Prevalence of Wolbachia')
            
            
            if p_m(end)>10^-2; lcolor = 'r'; else; lcolor = 'b'; end
            subplot(1,2,2); hold on            
            scatter(R0_wol*ones(size(p_m)),p_m,20,'MarkerEdgeColor',lcolor,'MarkerFaceColor',lcolor)
            hold off
            ylim([0 1])
            xlabel('R0 wol')
            ylabel('Prevalence of malaria')
            
            pause(0.1)
            
%             figure
%             plot(t,p_w,'linewidth',2)

        end
    end
end
toc



%% Plot output
% SH = y(:,1);
% EH = y(:,2);
% AH = y(:,3);
% DH = y(:,4);
% Ie = y(:,5);
% SU = y(:,6);
% EU = y(:,7);
% IU = y(:,8);
% SW = y(:,9);
% EW = y(:,10);
% IW = y(:,11);

% figure(1)
% subplot(2,2,1)
% plot(t,[SH EH AH DH],'linewidth',2)
% set(gca,'fontsize',18)
% xlabel('Time, days')
% ylabel('Population (human)')
% legend('S_H','E_H','A_H','D_H')
% ylim([0 SH0])
%
% subplot(2,2,2)
% plot(t,Ie,'linewidth',2)
% set(gca,'fontsize',18)
% xlabel('Time, days')
% ylabel('Immunity')
% legend('I_e')
%
% subplot(2,2,3)
% plot(t,[SU EU IU],'linewidth',2)
% set(gca,'fontsize',18)
% xlabel('Time, days')
% ylabel('Population (mosquito)')
% legend('S_U','E_U','I_U')
% ylim([0 P.Kf])
%
% subplot(2,2,4)
% plot(t,[SW EW IW],'linewidth',2)
% set(gca,'fontsize',18)
% xlabel('Time, days')
% ylabel('Population (mosquito)')
% legend('S_W','E_W','I_W')
% ylim([0 P.Kf])


% figure
% plot(t,[EH AH DH],'linewidth',2)
% set(gca,'fontsize',18)
% xlabel('Time, days')
% ylabel('Proportion of population')
% legend('E_H','A_H','D_H')

toc