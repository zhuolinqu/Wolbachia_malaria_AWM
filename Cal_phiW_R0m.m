function [r, R0m] = Cal_phiW_R0m(P)
% solve R0m(phiW)=1
options = optimset('TolX',10^-3);
r = fminbnd(@(x) (R0m_phiW(x,P)-1).^2,0,100,options);
R0m = R0m_phiW(r,P);
end

function R0m = R0m_phiW(x,P)

P.phiW = x;
SS_matW = EquilibriumState_w(P);
Wol_EEp = SS_matW(3,1:end-1);
Su = Wol_EEp(1); Sw = Wol_EEp(2);
R0m = Cal_R0_malaria(Su,Sw,P);

end