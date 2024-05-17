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
peak_list = NaN(size(t_malaria_control_grid));
rebound_list = NaN(size(t_malaria_control_grid));
dt = 1;
eff_malaria = 0.8; % efficacy = reduce malaria infection by X %
p = 0.5;
%% initial condition EE
SS_matM = EquilibriumState_m(P);
ySS = SS_matM(4,1:11);
options = odeset('AbsTol',1e-10,'RelTol',1e-10);
[tEE,yEE] = ode45(@BaseModel,[0 tinit],ySS,options,P);
yEE = yEE(end,:);
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
    % summarize output
    SH = y(:,1); EH = y(:,2); AH = y(:,3); DH = y(:,4); 
    NH = SH+EH+AH+DH; 
    p_malaria_inf_human = (AH+DH)./NH;
    [~,ind1] = min(abs(t-t_malaria_control));
    malaria_inf_start = p_malaria_inf_human(ind1);
    p_malaria_after = p_malaria_inf_human(ind1:end); 
    [~,ind2] = max(p_malaria_after);
    malaria_inf_max = p_malaria_after(ind2);
    peak_list(itime) = malaria_inf_max;
    rebound_list(itime) = malaria_inf_max-malaria_inf_start;
    t_mat(itime,:) = t; y_mat(itime,:,:) = y;
end
%% Plotting
figure_setups; hold on
plot(t_malaria_delay_grid,peak_list,'-','DisplayName','Peak of malaria since release')
plot(t_malaria_delay_grid,rebound_list,'--','DisplayName','Rebound of malaria')
xlabel('Time since \textit{Wolbachia} release, days')
ylabel('Fraction of infection $(A_H\,\&\,D_H)$')
legend

% figure_setups; 
% plot(t_mat(1,:),(y_mat(1,:,3)+y_mat(1,:,4))./sum(y_mat(1,:,1:4),3),'-')
