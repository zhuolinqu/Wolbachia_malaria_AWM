clearvars; 
% close all; 
clc
tic
%% Parameters
Baseline_params_malaria;
P = Baseline_params_stephensi(P);
P.vw = 0.95; P.vu = 1- P.vw;
%% Run model
% Time frame
tspan = [0 1000];

% Initial conditions
[R0w, G0w, G0u] = Cal_R0_wolbachia(P);

SS_mat = EquilibriumState_m(P);
% yinit(1:5) = SS_mat(4,1:5); % malaria endemic in human, no Wolbachia
% yinit(6:8) = SS_mat(4,6:8); % malaria endemic in mosquito, no Wolbachia
yinit = SS_mat(5,1:11);
% yinit(9) = yinit(9)*0.8; % SW
% yinit(10) = yinit(10)*0.8; % EW
yinit(11) = yinit(11)*0.8; % IW

% yinit = [SH0; EH0; AH0; DH0; Ie0; SU0; EU0; IU0; SW0; EW0; IW0];

options = odeset('AbsTol',1e-10,'RelTol',1e-10);
[t,y] = ode45(@BaseModel,tspan,yinit,options,P);

toc

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
NH = SH+EH+AH+DH;
NW = SW+EW+IW;
NU = SU+EU+IU;
NM = NW+NU;

%% Plotting
plot_solu_all_time; % plot solutions in time
% plot_EAD_only_time; % plot diseased groups only
% plot_sigmoids_time; % plot sigmoids
plot_infection_time; % plot wolbachia & malaria prevalence

