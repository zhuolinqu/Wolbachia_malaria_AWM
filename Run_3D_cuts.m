clearvars;
close all; clc
format long

%% Parameters
Baseline_params_malaria;
Baseline_params_stephensi;
P.vw = 1; P.vu = 1-P.vw;

phiW_list = linspace(0.5,5,10);

for iphi = 1:length(phiW_list)
    
    P.phiW = phiW_list(iphi);
    
    [R0w, G0w, G0u] = Cal_R0_wolbachia(P);
    a = P.vu/P.vw;
    b = (1-P.ci)/R0w+(P.vu/P.vw-1);
    c = (1-R0w)/R0w;
    
    if R0w>1
        case_id = 3; % always move to high stable EE
        rp = (-b + sqrt(b^2-4*a*c))/(2*a);
    elseif b^2-4*a*c>=0
        case_id = 2; % bistable
        if P.vw<1
            rn = (-b - sqrt(b^2-4*a*c))/(2*a);
            rp = (-b + sqrt(b^2-4*a*c))/(2*a);
        elseif P.vw==1
            rn = -c/b;
        end
    elseif b^2-4*a*c<0
        case_id = 1; % always move to DFE
    end
    
    tinit = 1000;
    tconti = 5000;
    %% Initial condition: Generate EE for malaria (no Wolbachia)
    SH0 = P.gH/P.muH-1; EH0 = 1; AH0 = 0; DH0 = 0;
    Ie0 = 0;
    SU0 = P.Kf; EU0 = 0; IU0 = 0; SW0 = 0; EW0 = 0; IW0 = 0;
    yinit = [SH0; EH0; AH0; DH0; Ie0; SU0; EU0; IU0; SW0; EW0; IW0];
    
    options = odeset('AbsTol',1e-10,'RelTol',1e-10);
    [t,y] = ode45(@BaseModel,linspace(0,tinit,50),yinit,options,P);
    
    yinit = y(end,:); % Wolbachia free malaria endemic state
    
    %% turn on Wolbachia
    %     if case_id==1 || case_id==3
    %         Nw_threshold = 1; % release one infected mosquito
    %         yinit(9) = Nw_threshold;
    %     elseif case_id==2
    %         r = rn*1.01;
    %         Nu_threshold = P.Kf/(1+r*P.mufu/P.mufw)*(1-1/G0w);
    %         Nw_threshold = r*P.mufu/P.mufw*Nu_threshold;
    %         % release above Wolbachia threshold
    %         temp = sum(yinit(6:8));
    %         yinit(6) = Nu_threshold*yinit(6)/temp; % Su
    %         yinit(7) = Nu_threshold*yinit(7)/temp; % Eu
    %         yinit(8) = Nu_threshold*yinit(8)/temp; % Iu
    %         yinit(9) = Nw_threshold;
    %     end
    
    
    frac_0 = 0.5; % initial infection fraction;
    Nw_0 = P.Kf*(1-1/G0w)*frac_0;
    Nu_0 = P.Kf*(1-1/G0w)*(1-frac_0);
    
    
    temp = sum(yinit(6:8));
    yinit(6) = Nu_0*yinit(6)/temp; % Su
    yinit(7) = Nu_0*yinit(7)/temp; % Eu
    yinit(8) = Nu_0*yinit(8)/temp; % Iu
    yinit(9) = Nw_0;
    
    [t1,y1] = ode45(@BaseModel,linspace(0,tconti,500),yinit,options,P);
    
    t = t1;
    y = y1;
    
    
%     t = [t(1:end-1);t1];
%     y = [y(1:end-1,:);y1];
    
    %% remove the redundant time steps after SS
    y_der = vecnorm(diff(y),2,2);
    ind = length(t) - find(flip(y_der)>700,1); % solution index when norm(change) < 500
    t = t(1:end-1);
    y = y(1:end-1,:);
    t(ind:end)=[];
    y(ind:end,:)=[];
    %% Plot output
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
    
    
    figure
    hold on
    col = t;
    scatter(NW./NM,(AH+DH)./NH,[],col,'o','filled')
    set(gca,'fontsize',18)
    xlabel('Wolbachia prevalence')
    ylabel('Malaria prevalence in humans')
    axis([0 1 0 1])
    colorbar
    caxis([0 300])
    title(['R0w=',num2str(R0w), '  case = ', num2str(case_id)])
    
%     if case_id==1
%         pW_final = 0;
%     elseif case_id==2 || case_id==3
%         if P.vw<1 & 
%             r = rp;
%             pW_final = (r*P.mufu/P.mufw)/(1+r*P.mufu/P.mufw);
%         elseif P.vw==1
%             pW_final=1;
%         end
%     end
    pW_final = NW(end)/NM(end);
    mW_final = (AH(end)+DH(end))/NH(end);
    plot(pW_final,mW_final,'*','linewidth',2)
    
end

% figure
% subplot(2,2,1)
% plot(t,[SH EH AH DH],'linewidth',2)
% set(gca,'fontsize',18)
% xlabel('Time, days')
% ylabel('Population (human)')
% legend('S_H','E_H','A_H','D_H')
% ylim([0 P.gH/P.muH])
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
%%
% figure
% plot(NW./NM,(AH+DH)./NH,'o')
% set(gca,'fontsize',18)
% xlabel('Wolbachia prevalence')
% ylabel('Malaria prevalence in humans')
%%
% figure
% plot(t,NM,'linewidth',2)
% set(gca,'fontsize',18)
% xlabel('Wolbachia prevalence')
% ylabel('Malaria prevalence in humans')

%%
% figure
% plot(t,[EH AH DH],'linewidth',2)
% set(gca,'fontsize',18)
% xlabel('Time, days')
% ylabel('Proportion of population')
% legend('E_H','A_H','D_H')
%%
% figure
% subplot(1,2,1)
% hold on
% plot(t,NU,'linewidth',2)
% plot(t,NU_diag,'linewidth',2)
% set(gca,'fontsize',18)
% xlabel('Time, days')
% ylabel('Population (mosquito)')
% legend('NU','NU diag')
%
% subplot(1,2,2)
% hold on
% plot(t,NW,'linewidth',2)
% plot(t,NW_diag,'linewidth',2)
% set(gca,'fontsize',18)
% xlabel('Time, days')
% ylabel('Population (mosquito)')
% legend('NW','NW diag')
% %%
% figure
% % subplot(1,2,1)
% hold on
% plot(t,(AH+DH)./NH,'linewidth',2)
% plot(t,NW./NM,'linewidth',2)
% set(gca,'fontsize',18)
% xlabel('Time, days')
% ylabel('Proportion of population')
% legend('Malaria prevalence in humans','Wolbachia prevalence')

% subplot(1,2,2)
% hold on
% plot(t,(AH+DH)./NH,'linewidth',2)
% plot(t,NW./NM,'linewidth',2)
% plot(t,NW_diag./(NW_diag+NU_diag),'linewidth',2)
% set(gca,'fontsize',18)
% xlabel('Time, days')
% ylabel('Proportion of population')
% legend('Malaria prevalence in humans','Wolbachia prevalence','diagnostic')

