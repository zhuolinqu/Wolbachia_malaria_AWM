function [R0, RHM, RMH] = Cal_R0_malaria(NH,NM,Sw,Su,P)
% basic reproduction number for malaria transmission, given the Wolbachia
% level Su and Sw.



BM = P.bm*P.bh*NH/(P.bm*NM+P.bh*NH);
BH = P.bm*P.bh*NM/(P.bm*NM+P.bh*NH);

% mosquito to human
RMH = BH*P.betaM*(P.alpha*Sw/(Su+Sw)*P.sigma/(P.mufw+P.sigma)/P.mufw+...
    Su/(Su+Sw)*P.sigma/(P.mufu+P.sigma)/P.mufu);

% human to mosquito
RHM = BM*P.h/(P.h+P.muH)*(P.rho*P.betaD/(P.rD+P.muD+P.muH)+...
                         (1-P.rho)*P.betaA/(P.rA+P.muH)+...
                         P.rho*(1-P.phi)*P.rD/(P.rD+P.muD+P.muH)*P.betaA/(P.rA+P.muH));

R0 = sqrt(RHM*RMH);

end