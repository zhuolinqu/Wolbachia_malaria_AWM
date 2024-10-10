%% Numerical simulations on mosquito control pre wolbachia release
clear all;
close all;
clc

%% Parameters
Baseline_params_malaria;
P = Baseline_params_stephensi(P);
t_Wolbachia_release = 150; % time of Wolbachia release 
t_mosquito_delay_grid = [-100,-50,-30]; % delay of mosquito control following W release
% t_mosquito_delay_grid = [-2,-1,0]; % delay of mosquito control following W release
t_mosquito_control_grid = t_Wolbachia_release+t_mosquito_delay_grid;
t_final = 600; % final time for simulation
total_inf_days_list = NaN(size(t_mosquito_control_grid));
eff_mosquito = 0.8; % efficacy = reduce mosquitoes by X %
p = 0.5;
options = odeset('AbsTol',1e-10,'RelTol',1e-10); dt = 0.5;
%% initial condition EE
SS_matM = EquilibriumState_m(P);
yEE = SS_matM(4,1:11);
%% Baseline scenario
tspan = 0:dt:t_Wolbachia_release;
[tbase,ybase] = ode45(@BaseModel,tspan,yEE,options,P);
y0 = Wol_release(ybase(end,:),p,yEE);
tspan = t_Wolbachia_release:dt:t_final;
[tbase2,ybase2] = ode45(@BaseModel,tspan,y0,options,P);
ybase = [ybase(1:end-1,:);ybase2]; tbase = [tbase(1:end-1,:);tbase2];
AH = ybase(:,3); DH = ybase(:,4);
total_inf_days_0 = Cal_inf_days(tbase,AH,DH);
%% for different mosquito control timing
for itime = 1:length(t_mosquito_control_grid)
    disp([num2str(itime),'/',num2str(length(t_mosquito_control_grid))])
    t_mosquito_control = t_mosquito_control_grid(itime);
    delay = t_mosquito_delay_grid(itime);
    if delay < 0 % mosquito control before W release
        tspan = 0:dt:t_mosquito_control;
        [t1,y1] = ode45(@BaseModel,tspan,yEE,options,P);
        t = t1; y = y1;
        y0 = mosquito_control(y1(end,:),eff_mosquito);
        tspan = t_mosquito_control:dt:t_Wolbachia_release;
        [t2,y2] = ode45(@BaseModel,tspan,y0,options,P);
        t = [t(1:end-1);t2]; y = [y(1:end-1,:);y2];
        y0 = Wol_release(y2(end,:),p,yEE);
        tspan = t_Wolbachia_release:dt:t_final;
        [tfinal,yfinal] = ode45(@BaseModel,tspan,y0,options,P);
        t = [t(1:end-1);tfinal];
        y = [y(1:end-1,:);yfinal];
    elseif delay == 0
        tspan = 0:dt:t_Wolbachia_release;
        [t1,y1] = ode45(@BaseModel,tspan,yEE,options,P);
        t = t1; y = y1;
        y0 = Wol_release(y1(end,:),p,yEE);
        y0 = mosquito_control(y0,eff_mosquito);
        tspan = t_Wolbachia_release:dt:t_final;
        [t2,y2] = ode45(@BaseModel,tspan,y0,options,P);
        t = [t(1:end-1,:);t2];
        y = [y(1:end-1,:);y2];
    end
    AH = y(:,3); DH = y(:,4);
    total_inf_days_list(itime) = Cal_inf_days(t,AH,DH);
    t_mat(itime,:) = t; y_mat(itime,:,:) = y;
end
%% Plotting
% figure_setups; hold on
% plot(t_malaria_delay_grid,total_inf_days_0-total_inf_days_list,'-o','DisplayName','\# of infection days reduced')
% xlabel('Time since Wolbachia release, days')
% ylabel('infected days reduced')
%% malaria human prevalence in time
colour_mat = [0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560];
h = figure_setups; 
yyaxis left
axis([0, 600, 0, 15000])
hold on
p_0 = plot(tbase,ybase(:,3)+ybase(:,4),':','Color',colour_mat(1,:),'DisplayName','no control');
label_linestyle = {'-','-.','--'};
p_1 = plot(t_mat(1,:),y_mat(1,:,3)+y_mat(1,:,4),'LineStyle',label_linestyle{1},'Color',colour_mat(2,:),'DisplayName',['delay = ',num2str(t_mosquito_delay_grid(1))]);
p_2 = plot(t_mat(2,:),y_mat(2,:,3)+y_mat(2,:,4),'LineStyle',label_linestyle{2},'Color',colour_mat(3,:),'DisplayName',['delay = ',num2str(t_mosquito_delay_grid(2))]);
p_3 = plot(t_mat(3,:),y_mat(3,:,3)+y_mat(3,:,4),'LineStyle',label_linestyle{3},'Color',colour_mat(4,:),'DisplayName',['delay = ',num2str(t_mosquito_delay_grid(3))]);
xlabel('time, days')
ylabel('Infectious population ($A_H+D_H$)')
ll = legend('Position',[0.61,0.61,0.285,0.25]); % delay = [-100,-50,-30];
% ll = legend('Position',[0.645,0.61,0.25,0.25]); % delay = [-2,-1,0];
getLegendPosition;
plotLegendBoundary(h, x0, x1, y0, y1, 'FaceColor', 'w');
fac = 0.03;
plotLegendLine(h, x0+dx*fac, y1-0.04 - 0.08*0, p_0, [1,0],  1);
plotLegendLine(h, x0+dx*fac, y1-0.04 - 0.08*1, p_1, [1,0],  1);
plotLegendLine(h, x0+dx*fac, y1-0.04 - 0.08*2, p_2, [1.0, 0.2, 0.5, 0.2],  9); % '-.'
plotLegendLine(h, x0+dx*fac, y1-0.04 - 0.08*3, p_3, [0.5, 0.2, 0.5,0.2], 8); %'--'
yyaxis right
axis([0, 600, 0, 1])
hold on
plot(tbase,sum(ybase(:,9:11),2)./sum(ybase(:,6:11),2),':','Color',colour_mat(1,:),'DisplayName','Wolbachia prev. in mosquitoes')
plot(t_mat(1,:),sum(y_mat(1,:,9:11),3)./sum(y_mat(1,:,6:11),3),'LineStyle',label_linestyle{1},'Color',colour_mat(2,:),'DisplayName','Wolbachia prev. in mosquitoes')
plot(t_mat(2,:),sum(y_mat(2,:,9:11),3)./sum(y_mat(2,:,6:11),3),'LineStyle',label_linestyle{2},'Color',colour_mat(3,:),'DisplayName','Wolbachia prev. in mosquitoes')
plot(t_mat(3,:),sum(y_mat(3,:,9:11),3)./sum(y_mat(3,:,6:11),3),'LineStyle',label_linestyle{3},'Color',colour_mat(4,:),'DisplayName','Wolbachia prev. in mosquitoes')
ylabel('\textit{Wolbachia} prevalence')
title(['Mosquito control efficacy = ',num2str(eff_mosquito)])