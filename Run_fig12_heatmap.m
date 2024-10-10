%% Numerical simulations on mosquito control pre wolbachia release
clear all;
% close all;
clc

%% Parameters
Baseline_params_malaria;
P = Baseline_params_stephensi(P);

t_mosquito_delay_grid = [-30:1:0]; % delay of mosquito control following W release
eff_mosquito_grid = [0.2:0.1:0.8]; % efficacy = reduce mosquitoes by X %

total_inf_days_list = NaN(length(t_mosquito_delay_grid),length(eff_mosquito_grid));
peak_list = NaN(length(t_mosquito_delay_grid),length(eff_mosquito_grid));
rebound_list = NaN(length(t_mosquito_delay_grid),length(eff_mosquito_grid));

t_Wolbachia_release = 150; % time of Wolbachia release
t_mosquito_control_grid = t_Wolbachia_release+t_mosquito_delay_grid;
t_final = 600; % final time for simulation
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
%% for different mosquito control timing and efficacy
for itime = 1:length(t_mosquito_control_grid)
    for jeff = 1:length(eff_mosquito_grid)
        disp([num2str(itime),'_',num2str(jeff),'/',num2str(length(t_mosquito_control_grid)),'_',num2str(length(eff_mosquito_grid))])
        t_mosquito_control = t_mosquito_control_grid(itime);
        delay = t_mosquito_delay_grid(itime);
        eff_mosquito = eff_mosquito_grid(jeff);
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
        % summarize output
        SH = y(:,1); EH = y(:,2); AH = y(:,3); DH = y(:,4);
        NH = SH+EH+AH+DH;
        p_malaria_inf_human = (AH+DH)./NH;
        [~,ind1] = min(abs(t-t_mosquito_control));
        malaria_inf_start = p_malaria_inf_human(ind1);
        p_malaria_after = p_malaria_inf_human(ind1+1:end);
        [~,ind2] = max(p_malaria_after);
        malaria_inf_max = p_malaria_after(ind2);
        peak_list(itime,jeff) = malaria_inf_max;
        rebound_list(itime,jeff) = malaria_inf_max-malaria_inf_start;
        AH = y(:,3); DH = y(:,4);
        total_inf_days_list(itime,jeff) = Cal_inf_days(t,AH,DH);

        % t_mat(itime,jeff,:) = t; y_mat(itime,jeff,:,:) = y;
    end
end
%% Plot (control time, control efficacy, QOI)
QOI = total_inf_days_list;
figure_setups; hold on; grid off
days_reduced = total_inf_days_0-total_inf_days_list;
zz = days_reduced;
[xx,yy] = meshgrid(t_mosquito_delay_grid,eff_mosquito_grid);
[C,h] = contourf(xx,yy,log10(zz'),25,'edgecolor','none');
shading interp
xlabel({'Deployment time of mosquito control, days','(relative to \textit{Wolbachia} release)'})
ylabel('Mosquito control efficacy')
title('Infectious days reduced');
set(gca,'YDir','normal');
colormap jet
colorbar
clim([4 6])
Ticks = [4:6];
TickLabels = {'10^{4}','10^{5}','10^{6}'};
colorbar('Ticks',Ticks,'TickLabels',TickLabels)
xlim([min(t_mosquito_delay_grid),max(t_mosquito_delay_grid)])
ylim([min(eff_mosquito_grid),max(eff_mosquito_grid)])

%% Infectious days reduced vs. deployment time @ efficacy = 0.8;
figure_setups; hold on
plot(t_mosquito_delay_grid,zz(:,end));
xlabel({'Deployment time of mosquito control, days','(relative to \textit{Wolbachia} release)'})
ylabel('Infectious days reduced')
ylim([0 12*10^5])
grid on
plot(t_mosquito_delay_grid(end),zz(end,end),'md')
title(['~~~~~~Mosquito control efficacy = ',num2str(eff_mosquito_grid(end))])
