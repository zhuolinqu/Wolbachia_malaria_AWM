clearvars; 
% close all; 
clc

tic
%% Parameters
Baseline_params_malaria;
P = Baseline_params_stephensi(P);
P.phiW = 0.7;
P.vw = 0.95; P.vu = 1- P.vw;
%% Run model
% Time frame
tspan = [0 1000];

% Initial conditions
[R0w, G0w, G0u] = Cal_R0_wolbachia(P);

SS_mat = EquilibriumState_m(P);
yinit = SS_mat(6,1:11);
% SH0 = 0.5*P.gH/P.muH;
% EH0 = 0.5*P.gH/P.muH;
% AH0 = 0;
% DH0 = 0;
% Ie0 = 0;
% SU0 = 1;
% EU0 = 0;
% IU0 = 0;
% SW0 = 0;% P.Kf*(1-1/G0w); 
% EW0 = 0;
% IW0 = 0;
% yinit(9) = yinit(9)*0.8;
% yinit(10) = yinit(10)*0.8;
% yinit(11) = yinit(11)*0.8;

% yinit = [SH0; EH0; AH0; DH0; Ie0; ...
%         SU0; EU0; IU0; ...
%         SW0; EW0; IW0];

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

%% plot solutions in time
figure_setups;
subplot(2,2,1)
plot(t,[SH EH AH DH])
xlabel('Time, days')
ylabel('Population (human)')
legend('$S_H$','$E_H$','$A_H$','$D_H$')
ylim([0 P.gH/P.muH])

subplot(2,2,2)
plot(t,Ie./NH)
xlabel('Time, days')
ylabel('Immunity')
legend('$I_e$')

subplot(2,2,3)
plot(t,[SU EU IU])
xlabel('Time, days')
ylabel('Population (mosquito)')
legend('$S_U$','$E_U$','$I_U$')
ylim([0 P.Kf])

subplot(2,2,4)
plot(t,[SW EW IW])
xlabel('Time, days')
ylabel('Population (mosquito)')
legend('$S_W$','$E_W$','$I_W$')
ylim([0 P.Kf])

%% plot diseased groups only
figure_setups;
subplot(1,2,1)
plot(t,[EH AH DH])
xlabel('Time, days')
ylabel('Population (human)')
legend('$E_H$','$A_H$','$D_H$')
subplot(1,2,2)
plot(t,[EH AH DH]./NH)
xlabel('Time, days')
ylabel('Proportion (human)')
legend('$E_H$','$A_H$','$D_H$')

%% plot wolbachia prevalence
figure_setups;
plot(t,NW./(NW+NU))
xlabel('Time, days')
ylabel('Wolbachia prevalence')
%% plot sigmoids
figure_setups;
[rho, phi, psi] = sigmoid_prob(Ie./NH, P);
plot(t,rho,t,phi,t,psi)
legend('rho','phi','psi')