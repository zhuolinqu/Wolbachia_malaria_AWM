function phiW = Cal_phiW_R0w(P)
% solve R0w(phiW) = 1

phiW = (P.mufw*P.phiU)/(P.mufu*P.vw);
