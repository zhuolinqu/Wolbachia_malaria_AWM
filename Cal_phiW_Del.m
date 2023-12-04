function phiW = Cal_phiW_Del(P)
% solve Delta(phiW) = 0
% Delta = b^2 - 4ac with
%   a = vu/vw;
%   b = (1 - ci)/R0w + (vu/vw - 1);
%   c = (1 - R0w)/R0w;

undersqrt = sqrt(P.ci*P.mufu^2*P.mufw^2*P.phiU^2*(-1+P.vw)*(-1+P.ci*P.vw));

phiW1 = (P.mufu*P.mufw*P.phiU*(1+P.ci-2*P.ci*P.vw)-2*undersqrt)/P.mufu^2;
phiW2 = (P.mufu*P.mufw*P.phiU*(1+P.ci-2*P.ci*P.vw)+2*undersqrt)/P.mufu^2;

phiW = max(phiW1,phiW2);
