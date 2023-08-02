function dy = BaseModel_debug(t,y,P)

SH = y(1);
EH = y(2);
AH = y(3);
DH = y(4);
Ie = y(5);
SU = y(6);
EU = y(7);
IU = y(8);
SW = y(9);
EW = y(10);
IW = y(11);
NU_diag = y(12);
NW_diag = y(13);

NH = SH + EH + AH + DH;
NU = SU + EU + IU;
NW = SW + EW + IW;
NM = NU + NW;

BM = P.bm*P.bh*NH/(P.bm*NM+P.bh*NH);
BH = P.bm*P.bh*NM/(P.bm*NM+P.bh*NH);

LamH = BH*P.betaM*(IU+IW)/NM;
LamM = BM*(P.betaD*DH+P.betaA*AH)/NH;

psi = P.psi;
rho = P.rho;
phi = P.phi;

f_LamH = LamH/(P.gamma*LamH+1);

dSH = P.gH - LamH*SH + phi*P.rD*DH ...
    + P.rA*AH - P.muH*SH;
dEH = LamH*SH - P.h*EH - P.muH*EH;
dAH = rho*P.h*EH - (1-psi)*LamH*AH ...
    + (1-phi)*P.rD*DH - P.rA*AH - P.muH*AH;
dDH = (1-rho)*P.h*EH + (1-psi)*LamH*AH - P.rD*DH ...
    - (P.muH+P.muD)*DH;
dIe = f_LamH*(P.cS*SH+P.cE*EH+P.cA*AH+P.cD*DH) ...
    - (1/P.de + P.muH + P.muD*DH/NH)*Ie;

gU = P.bf*P.phiU*NU/(NU+P.mufw/P.mufu*NW)*(1-(NU+NW)/P.Kf)*NU...
    + (1-P.ci)*P.bf*P.phiU*P.mufw/P.mufu*NW/(NU+P.mufw/P.mufu*NW)*(1-(NU+NW)/P.Kf)*NU...
    + P.vu*P.bf*P.phiW*(1-(NU+NW)/P.Kf)*NW;
gW = P.vw*P.bf*P.phiW*(1-(NU+NW)/P.Kf)*NW;

dSU = -LamM*SU + gU - P.mufu*SU;
dEU = LamM*SU - P.sigma*EU - P.mufu*EU;
dIU = P.sigma*EU - P.mufu*IU;
dSW = -P.alpha*LamM*SW + gW - P.mufw*SW;
dEW = P.alpha*LamM*SW - P.sigma*EW - P.mufw*EW;
dIW = P.sigma*EW - P.mufw*IW;

gU_diag = P.bf*P.phiU*NU_diag/(NU_diag+P.mufw/P.mufu*NW_diag)*(1-(NU_diag+NW_diag)/P.Kf)*NU_diag...
    + (1-P.ci)*P.bf*P.phiU*P.mufw/P.mufu*NW_diag/(NU_diag+P.mufw/P.mufu*NW_diag)*(1-(NU_diag+NW_diag)/P.Kf)*NU_diag...
    + P.vu*P.bf*P.phiW*(1-(NU_diag+NW_diag)/P.Kf)*NW_diag;
gW_diag = P.vw*P.bf*P.phiW*(1-(NU_diag+NW_diag)/P.Kf)*NW_diag;

dNU_diag =  gU_diag - P.mufu*NU_diag;
dNW_diag =  gW_diag- P.mufw*NW_diag;

dy = [dSH; dEH; dAH; dDH; dIe; dSU; dEU; dIU; dSW; dEW; dIW; dNU_diag; dNW_diag];

end