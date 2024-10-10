%% Numerical simulations on malaria control pre or post w/ wolbachia release
clearvars;
close all;
clc

%% Parameters
Baseline_params_malaria;
P = Baseline_params_stephensi(P);
tinit = 10; % start of simulation
t_Wolbachia_release = 150; % time of Wolbachia release 
t_malaria_delay_grid = [-50, 40, 100]; % delay of malaria control following W release
t_malaria_control_grid = t_Wolbachia_release+t_malaria_delay_grid;
t_final = 600; % final time for simulation
total_inf_days_list = NaN(size(t_malaria_control_grid));
dt = 1;
eff_malaria = 0.8; % efficacy = treat malaria infection by X %
p = 0.5;
%% initial condition EE
SS_matM = EquilibriumState_m(P);
ySS = SS_matM(4,1:11);
options = odeset('AbsTol',1e-10,'RelTol',1e-10);
[tEE,yEE] = ode45(@BaseModel,[0 tinit],ySS,options,P);
yEE = yEE(end,:);
%% Baseline scenario
tspan = tinit:dt:t_Wolbachia_release;
[tbase,ybase] = ode45(@BaseModel,tspan,yEE,options,P);
y0 = Wol_release(ybase(end,:),p);
tspan = t_Wolbachia_release:dt:t_final;
[tbase2,ybase2] = ode45(@BaseModel,tspan,y0,options,P);
ybase = [ybase(1:end-1,:);ybase2];
tbase = [tbase(1:end-1,:);tbase2];
AH = ybase(:,3); DH = ybase(:,4);
total_inf_days_0 = Cal_inf_days(tbase,AH,DH);
%% for different malaria control timing
for itime = 1:length(t_malaria_control_grid)
    disp([num2str(itime),'/',num2str(length(t_malaria_control_grid))])
    t_malaria_control = t_malaria_control_grid(itime);
    delay = t_malaria_delay_grid(itime);
    if delay < 0 % malaria control before W release
        tspan = tinit:dt:t_malaria_control;
        [t1,y1] = ode45(@BaseModel,tspan,yEE,options,P);
        t = t1; y = y1;
        y0 = malaria_control(y1(end,:),eff_malaria);
        tspan = t_malaria_control:dt:t_Wolbachia_release;
        [t2,y2] = ode45(@BaseModel,tspan,y0,options,P);
        t = [t(1:end-1);t2]; y = [y(1:end-1,:);y2];
        y0 = Wol_release(y2(end,:),p);
        tspan = t_Wolbachia_release:dt:t_final;
        [tfinal,yfinal] = ode45(@BaseModel,tspan,y0,options,P);
        t = [t(1:end-1);tfinal];
        y = [y(1:end-1,:);yfinal];
    elseif delay>0 % malaria control after W release
        tspan = tinit:dt:t_Wolbachia_release;
        [t1,y1] = ode45(@BaseModel,tspan,yEE,options,P);
        t = t1; y = y1;
        y0 = Wol_release(y1(end,:),p);
        tspan = t_Wolbachia_release:dt:t_malaria_control;
        [t2,y2] = ode45(@BaseModel,tspan,y0,options,P);
        t = [t(1:end-1);t2];  y = [y(1:end-1,:);y2];
        y0 = malaria_control(y2(end,:),eff_malaria);
        tspan = t_malaria_control:dt:t_final;
        [tfinal,yfinal] = ode45(@BaseModel,tspan,y0,options,P);
        t = [t(1:end-1);tfinal];
        y = [y(1:end-1,:);yfinal];
    elseif delay == 0 
        tspan = tinit:dt:t_Wolbachia_release;
        [t1,y1] = ode45(@BaseModel,tspan,yEE,options,P);
        t = t1; y = y1;
        y0 = Wol_release(y1(end,:),p);
        y0 = malaria_control(y0,eff_malaria);
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
%%
colour_mat = [0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560];
figure_setups; 
yyaxis left
hold on
plot(tbase,ybase(:,3)+ybase(:,4),'LineStyle',':','Color',colour_mat(1,:),'DisplayName','no control')
label_linestyle = {'-.','--'};
plot(t_mat(1,:),y_mat(1,:,3)+y_mat(1,:,4),'LineStyle',label_linestyle{1},'Color',colour_mat(3,:),'DisplayName',['delay = ',num2str(t_malaria_delay_grid(1))])
plot(t_mat(2,:),y_mat(2,:,3)+y_mat(2,:,4),'m-','DisplayName',['delay = ',num2str(t_malaria_delay_grid(2))])
plot(t_mat(3,:),y_mat(3,:,3)+y_mat(3,:,4),'LineStyle',label_linestyle{2},'Color',colour_mat(2,:),'DisplayName',['delay = ',num2str(t_malaria_delay_grid(3))])
xlabel('time, days')
ylabel('Infectious population ($A_H+D_H$)')
legend('Position',[0.59,0.61,0.26555602169038,0.25])
yyaxis right
plot(tbase,sum(ybase(:,9:11),2)./sum(ybase(:,6:11),2),'k-','DisplayName','Wolbachia prev. in mosquitoes')
set(gca,'YColor','k')
ylabel('\textit{Wolbachia} prevalence')






