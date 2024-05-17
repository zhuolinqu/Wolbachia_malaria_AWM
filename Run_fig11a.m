clearvars;
close all;
clc

%% Parameters
Baseline_params_malaria;
P = Baseline_params_stephensi(P);
tinit = 10; % start of simulation
t_Wolbachia_release = 150; % time of Wolbachia release 
t_malaria_delay_grid = [-50:2:100]; % delay of malaria control following W release
t_malaria_control_grid = t_Wolbachia_release+t_malaria_delay_grid;
t_final = 600; % final time for simulation
total_inf_days_list = NaN(size(t_malaria_control_grid));
dt = 1;
eff_malaria = 0.8; % efficacy = reduce malaria infection by X %
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
figure_setups; hold on
days_reduced = total_inf_days_0-total_inf_days_list;
plot(t_malaria_delay_grid,days_reduced,'-','DisplayName','\# of infection days reduced')
[~,ind] = max(days_reduced);
plot(t_malaria_delay_grid(ind),days_reduced(ind),'md')
xlabel('Time since \textit{Wolbachia} release, days')
ylabel('Infectious days reduced')