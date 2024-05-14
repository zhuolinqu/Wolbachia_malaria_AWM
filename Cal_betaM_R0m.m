function [betaM, R0m] = Cal_betaM_R0m(P,R0m_star)
% solve R0m(betaM)=R0m_star
options = optimset('TolX',10^-10);
betaM = fminbnd(@(x) (R0m_betaM(x,P)-R0m_star).^2,0,500,options);
R0m = R0m_betaM(betaM,P);
if (R0m-R0m_star)^2>1e-8
    % fail to find betaM root
    betaM = NaN;
    R0m = NaN;
end
end

function R0m = R0m_betaM(x,P)

P.betaM = x;
SS_matW = EquilibriumState_w(P);
Wol_DFE = SS_matW(1,1:end-1);
Su = Wol_DFE(1); Sw = Wol_DFE(2);
R0m = Cal_R0_malaria(Su,Sw,P);

end