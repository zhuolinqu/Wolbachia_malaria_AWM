clearvars; 
% close all; 
clc

%% Parameters
Baseline_params_malaria;
P = Baseline_params_stephensi(P);
P.vw = 0.95; P.vu = 1- P.vw;
tinit = 100; % [start, pre-release control]
tconti_pre = 20; % [pre-release control, W release]
tconti_release = 100; % [W release, post-release control]
tconti_post = 1000-tinit-tconti_pre-tconti_release; % [post-release control, end]
flag_malaria_control_pre = 0; % malaria control as pre-release
flag_malaria_control_post = 1; % malaria control as post-release
eff_mosquito = 0.8; % 0.8 efficacy = reduce wild mosquito pop by X %
p = 0.5;
eff_malaria = 0.8; % 0.8 efficacy = reduce malaria infection by X %
%% Baseline scenario
tspan = 0:tinit;
SS_matM = EquilibriumState_m(P);
y0 = SS_matM(4,1:11);
options = odeset('AbsTol',1e-10,'RelTol',1e-10);
[t,y] = ode45(@BaseModel,tspan,y0,options,P);
%% Pre-release intervention 
tspan = t(end):t(end)+tconti_pre;
y0 = y(end,:);
y0 = mosquito_control(y0,eff_mosquito);
if flag_malaria_control_pre==1
    y0 = malaria_control(y0,eff_malaria);
end
if tconti_pre>0
    [tt,yy] = ode45(@BaseModel,tspan,y0,options,P);
    t = [t(1:end-1);tt];
    y = [y(1:end-1,:);yy];
    y0 = y(end,:);
end
%% Wolbachia release 
tspan = t(end):t(end)+tconti_release;
% SS_matW = EquilibriumState_w(P);
% thres is around 49.82% = SS_matW(2,2)/(SS_matW(2,1)+SS_matW(2,2))*0.95; % Wolbachia threshold level
y0 = Wol_release(y0,p);
[tt,yy] = ode45(@BaseModel,tspan,y0,options,P);
t = [t(1:end-1);tt];
y = [y(1:end-1,:);yy];
%% post-release control
tspan = t(end):t(end)+tconti_post;
y0 = y(end,:);
if flag_malaria_control_post==1
    y0 = malaria_control(y0,eff_malaria);
end
[tt,yy] = ode45(@BaseModel,tspan,y0,options,P);
t = [t(1:end-1);tt];
y = [y(1:end-1,:);yy];

% summarize output
SH = y(:,1); EH = y(:,2); AH = y(:,3); DH = y(:,4); Ie = y(:,5);
SU = y(:,6); EU = y(:,7); IU = y(:,8); 
SW = y(:,9); EW = y(:,10); IW = y(:,11);
NH = SH+EH+AH+DH; NW = SW+EW+IW; NU = SU+EU+IU; NM = NW+NU;

%% Plotting
% plot_solu_all_time; % plot solutions in time
% plot_EAD_only_time; % plot diseased groups only
% plot_sigmoids_time; % plot sigmoids
plot_infection_time; % plot wolbachia & malaria prevalence

% saveas(gcf,['Results/release_malaria_pre1.eps'],'epsc')
% saveas(gcf,['Results/release_malaria_post.eps'],'epsc')
