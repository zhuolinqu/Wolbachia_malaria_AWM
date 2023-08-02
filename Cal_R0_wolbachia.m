function [R0, G0w, G0u] = Cal_R0_wolbachia(P)
% basic reproduction number for Wolbachia transmission in mosquito


R0 = P.vw*P.mufu*P.phiW/(P.mufw*P.phiU);
G0w = P.vw*P.bf*P.phiW/P.mufw;
G0u = P.bf*P.phiU/P.mufu;

end