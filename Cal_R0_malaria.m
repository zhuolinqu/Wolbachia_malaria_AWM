function [R0, RHM, RMH] = Cal_R0_malaria(Su,Sw,P)
% basic reproduction number for malaria transmission, given the Wolbachia
% level Su and Sw.

NH = P.gH/P.muH;
NM = Su+Sw;
BM = P.bm*P.bh*NH/(P.bm*NM+P.bh*NH);
BH = P.bm*P.bh*NM/(P.bm*NM+P.bh*NH);

[rho, phi, ~] = sigmoid_prob(0, P); % at malaria DFE, no exposure-acquired immunity

% mosquito to human
RMH = BH*P.betaM*(P.alpha*Sw/(Su+Sw)*P.sigma/(P.mufw+P.sigma)/P.mufw+...
    Su/(Su+Sw)*P.sigma/(P.mufu+P.sigma)/P.mufu);

% human to mosquito
RHM = BM*P.h/(P.h+P.muH)*((1-rho)*P.betaD/(P.rD+P.muD+P.muH)+...
                         rho*P.betaA/(P.rA+P.muH)+...
                         (1-rho)*(1-phi)*P.rD/(P.rD+P.muD+P.muH)*P.betaA/(P.rA+P.muH));

R0 = sqrt(RHM*RMH);

end