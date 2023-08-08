clearvars; close all; clc

tic
%% Parameters
Baseline_params_malaria;
P = Baseline_params_stephensi(P);

%% Run model
% Time frame
tspan = [0 500];

% Initial conditions
SH0 = 0;
EH0 = P.gH/P.muH;
AH0 = 0;
DH0 = 0;
Ie0 = 0;
SU0 = 0;
EU0 = 0;
IU0 = 0;
SW0 = P.Kf;
EW0 = 0;
IW0 = 0;

yinit = [SH0; EH0; AH0; DH0; Ie0; ...
        SU0; EU0; IU0; ...
        SW0; EW0; IW0];

options = [];
% options = odeset('NonNegative',1:11); % if we want to enforce non-negativity

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

figure(1)
subplot(2,2,1)
plot(t,[SH EH AH DH],'linewidth',2)
set(gca,'fontsize',18)
xlabel('Time, days')
ylabel('Population (human)')
legend('S_H','E_H','A_H','D_H')
ylim([0 P.gH/P.muH])

subplot(2,2,2)
plot(t,Ie,'linewidth',2)
set(gca,'fontsize',18)
xlabel('Time, days')
ylabel('Immunity')
legend('I_e')

subplot(2,2,3)
plot(t,[SU EU IU],'linewidth',2)
set(gca,'fontsize',18)
xlabel('Time, days')
ylabel('Population (mosquito)')
legend('S_U','E_U','I_U')
ylim([0 P.Kf])

subplot(2,2,4)
plot(t,[SW EW IW],'linewidth',2)
set(gca,'fontsize',18)
xlabel('Time, days')
ylabel('Population (mosquito)')
legend('S_W','E_W','I_W')
ylim([0 P.Kf])
%%
figure
plot(t,[EH AH DH],'linewidth',2)
set(gca,'fontsize',18)
xlabel('Time, days')
ylabel('Proportion of population')
legend('E_H','A_H','D_H')
