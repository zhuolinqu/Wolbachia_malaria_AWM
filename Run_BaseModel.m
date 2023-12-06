clearvars; 
% close all; 
clc

%% Parameters
Baseline_params_malaria;
P = Baseline_params_stephensi(P);
P.vw = 0.95; P.vu = 1- P.vw;
%% Initial condition
tspan = [0 1000];
SS_matM = EquilibriumState_m(P);
SS_matW = EquilibriumState_w(P);
thres = SS_matW(2,2)/(SS_matW(2,1)+SS_matW(2,2)); % Wolbachia threshold level
yinit = SS_matM(4,1:11);
p = thres*0.95;
% reduce malaria infection
% SH EH AH DH
yinit(1) = yinit(1)+yinit(2)/2+yinit(3)/2+yinit(4)/2;
yinit(2) = yinit(2)/2;
yinit(3) = yinit(3)/2;
yinit(4) = yinit(4)/2;
% SU EU IU
yinit(6) = yinit(6)+yinit(7)+yinit(8);
yinit(7) = yinit(7)-yinit(7);
yinit(8) = yinit(8)-yinit(8);
% SW EW IW
yinit(9) = (yinit(6)+yinit(7)+yinit(8))*p/(1-p); % SW --> release W-infected mosquitoes
yinit(10) = 0;
yinit(11) = 0; 
%% ODE Solver
options = odeset('AbsTol',1e-10,'RelTol',1e-10);
[t,y] = ode45(@BaseModel,tspan,yinit,options,P);

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

