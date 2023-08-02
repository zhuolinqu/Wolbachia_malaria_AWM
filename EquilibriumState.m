function [Sustar, Swstar] = EquilibriumState(r,P)

[R0w, G0w, G0u] = Cal_R0_wolbachia(P);
Sustar = P.Kf/(1+r*P.mufu/P.mufw)*(1-1/G0w);
Swstar = r*P.mufu/P.mufw*Sustar;
end