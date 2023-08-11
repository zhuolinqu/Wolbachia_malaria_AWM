clear all;close all;
% clc

Baseline_params_malaria;
P = Baseline_params_stephensi(P);
% P = Baseline_params_gambiae(P);
[R0w, G0w, G0u] = Cal_R0_wolbachia(P);

Sw = 0;
Su = P.Kf*(1-1/G0u);

[R0, RHM, RMH] = Cal_R0_malaria(Su,Sw,P);

R0