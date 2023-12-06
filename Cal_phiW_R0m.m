function [phiW, R0m] = Cal_phiW_R0m(P)
% solve R0m(phiW)=1
options = optimset('TolX',10^-10);
phiW = fminbnd(@(x) (R0m_phiW(x,P)-1).^2,0,500,options);
R0m = R0m_phiW(phiW,P);
if (R0m-1)^2>1e-8
    % R0m < 1 for all phiW
    phiW = NaN;
    R0m = NaN;
end
end

function R0m = R0m_phiW(x,P)

P.phiW = x;
SS_matW = EquilibriumState_w(P);
Wol_EEp = SS_matW(3,1:end-1);
Su = Wol_EEp(1); Sw = Wol_EEp(2);
R0m = Cal_R0_malaria(Su,Sw,P);

end