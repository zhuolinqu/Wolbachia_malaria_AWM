% Baseline_params_malaria
% P = Baseline_params_stephensi(P);
% P.vw = 0.99;
% P.betaM = P.betaM/2;
%[r, R0m] = Cal_phiW_R0m(P)
% Baseline_params_malaria
% P = Baseline_params_stephensi(P);
% P.phiW
P.vw = 0.85;
P.vu = 1-P.vw;

phiW_vec = 0:.1:500;
for j = 1:length(phiW_vec)
    x = phiW_vec(j);
R0m = R0m_phiW(x,P);
R(j) = R0m;
end
figure(10);plot(phiW_vec,R)

%%
Baseline_params_malaria
P = Baseline_params_stephensi(P);
P.phiW
P.vw = 0.91;
P.vu = 1-P.vw;
[phiW, R0m] = Cal_phiW_R0m(P)
P.phiW = phiW;
[R0wB, G0w, G0u] = Cal_R0_wolbachia(P);
R0wB
        

function R0m = R0m_phiW(x,P)

P.phiW = x;
SS_matW = EquilibriumState_w(P);
Wol_EEp = SS_matW(3,1:end-1);
Su = Wol_EEp(1); Sw = Wol_EEp(2);
R0m = Cal_R0_malaria(Su,Sw,P);

end