function  SS_mat = EquilibriumState_m(P,SS_mat_old)

% malaria-Wolbachia

% (row 1) DFE-DFE: no malaria and no Wolbachia
% (row 2) DFE-EE-: no malaria and unstable Wolbachia endemic
% (row 3) DFE-EE+: no malaria and stable Wolbachia endemic
% (row 4) EE-DFE: malaria endemic and no Wolbachia
% (row 5) EE-EE-: malaria endemic and unstable Wolbachia endemic
% (row 6) EE-EE+: malaria endemic and stable Wolbachia endemic

% SS_mat(:,1:11) solution (1:5) malaria + (6:11) wolbachia
% SS_mat(:,n=12) indicates stability for malaria. 1 = stable, 0 = unstable
% SS_mat(:,n=13) indicates stability for wolbachia. 1 = stable, 0 = unstable
% SS_mat(:,n=14) last column indicates stability overall. 1 = stable, 0 = unstable  

%[SH, EH, AH, DH, Ie, SU, EU, IU, SW, EW, IW]

SS_mat = NaN(6,14);
index_m = 12;
index_w = 13;
index_all = 14;

[R0w, G0w, G0u] = Cal_R0_wolbachia(P);

if exist('SS_mat_old','var')
    % if the nearby steady state is provided, use these values as initial
    % conditions for steady state that are numerically captured
    flag_nearby = 1;
else
    flag_nearby = 0;
end

%% Check
if G0u<1 || G0w<1
    disp('mosquito extinction')
    % mosquito extinction, check parameter range!
end

SS_matW = EquilibriumState_w(P);
Wol_DFE = SS_matW(1,1:end-1); Wol_DFE_sta = SS_matW(1,end);
Wol_EEm = SS_matW(2,1:end-1); Wol_EEm_sta = SS_matW(2,end);
Wol_EEp = SS_matW(3,1:end-1); Wol_EEp_sta = SS_matW(3,end);

%% (row 1) DFE-DFE: no malaria and no Wolbachia
SU = Wol_DFE(1); SW = Wol_DFE(2); EU = 0; IU = 0; EW = 0; IW = 0;
R0m = Cal_R0_malaria(SU,0,P);
% keyboard
SH = P.gH/P.muH; EH = 0; AH = 0; DH = 0; Ie = 0;
SS_mat(1,1:11) = [SH, EH, AH, DH, Ie, SU, EU, IU, SW, EW, IW];

if R0m<1
    SS_mat(1,index_m) = 1;
else
    SS_mat(1,index_m) = 0; 
end
SS_mat(1,index_w) = Wol_DFE_sta; 
SS_mat(1,index_all) = SS_mat(1,index_w)*SS_mat(1,index_m);

% if R0m<1 && Wol_DFE_sta==1
%     SS_mat(1,end) = 1;
% else
%     SS_mat(1,end) = 0;
% end

%% (row 2) DFE-EE-: no malaria and unstable Wolbachia endemic
SU = Wol_EEm(1); SW = Wol_EEm(2);
if ~isnan(SU) % if Wolbachia EE- exist, always unstable
    EU = 0; IU = 0; EW = 0; IW = 0;
    SH = P.gH/P.muH; EH = 0; AH = 0; DH = 0; Ie = 0;
    SS_mat(2,1:11) = [SH, EH, AH, DH, Ie, SU, EU, IU, SW, EW, IW];   
    R0m = Cal_R0_malaria(SU, SW, P);

    if R0m<1
        SS_mat(2,index_m) = 1;
    else
        SS_mat(2,index_m) = 0;
    end
    if Wol_EEm_sta~=0; keyboard;end % wolbachia EE- should always be unstable
    SS_mat(2,index_w) = Wol_EEm_sta;
    SS_mat(2,index_all) = SS_mat(2,index_w)*SS_mat(2,index_m);

    % SS_mat(2,end) = 0;
end

%% (row 3) DFE-EE+: no malaria and high Wolbachia endemic
SU = Wol_EEp(1); SW = Wol_EEp(2);
if ~isnan(SU) % if Wolbachia EE+ exist
    EU = 0; IU = 0; EW = 0; IW = 0;
    SH = P.gH/P.muH; EH = 0; AH = 0; DH = 0; Ie = 0;
    R0m = Cal_R0_malaria(SU, SW, P);
    SS_mat(3,1:11) = [SH, EH, AH, DH, Ie, SU, EU, IU, SW, EW, IW];
    
    if R0m<1
        SS_mat(3,index_m) = 1;
    else
        SS_mat(3,index_m) = 0;
    end
    SS_mat(3,index_w) = Wol_EEp_sta;
    SS_mat(3,index_all) = SS_mat(3,index_w)*SS_mat(3,index_m);

    % if R0m<1 && Wol_EEp_sta==1
    %     SS_mat(3,end) = 1;
    % else
    %     SS_mat(3,end) = 0;
    % end
end

%% (row 4) EE-DFE: malaria endemic and no Wolbachia
SU = Wol_DFE(1); SW = Wol_DFE(2); % SW = 0;
R0m = Cal_R0_malaria(SU,SW,P);

if R0m>1
    SH0 = P.gH/P.muH-1; EH0 = 1; AH0 = 0; DH0 = 0; Ie0 = 0;
    SU0 = SU; EU0 = 0; IU0 = 0; SW0 = 0; EW0 = 0; IW0 = 0;
    yinit = [SH0; EH0; AH0; DH0; Ie0; SU0; EU0; IU0; SW0; EW0; IW0];
    options = odeset('AbsTol',1e-10,'RelTol',1e-10);
    [~,y] = ode45(@BaseModel,linspace(0,10^4,50),yinit,options,P);
    SS_mat(4,1:11) = y(end,:);

    SU_frac = y(end,6)/sum(y(end,6:8),2);
    EU_frac = y(end,7)/sum(y(end,6:8),2);
    IU_frac = y(end,8)/sum(y(end,6:8),2);

    SS_mat(4,index_m) = 1;
    SS_mat(4,index_w) = Wol_DFE_sta;
    SS_mat(4,index_all) = SS_mat(4,index_w)*SS_mat(4,index_m);

    % if Wol_DFE_sta==1
    %     SS_mat(4,end) = 1;
    % else
    %     SS_mat(4,end) = 0;
    % end

end

%% (row 6) EE-EE+: malaria endemic and high Wolbachia endemic
NU = Wol_EEp(1); NW = Wol_EEp(2);
R0m = Cal_R0_malaria(NU,NW,P);
if R0m>1 && ~isnan(NU) % if Wolbachia EE+ exist
    if flag_nearby && ~isnan(sum(SS_mat_old(6,1:11),2))
        yinit = SS_mat_old(6,1:11);
    else
        SH0 = P.gH/P.muH*0; EH0 = P.gH/P.muH*0.3; AH0 = P.gH/P.muH*0.35; DH0 = P.gH/P.muH*0.35; Ie0 = 0;
        SU0 = NU*SU_frac; EU0 = NU*EU_frac; IU0 = NU*IU_frac;
        SW0 = NW*SU_frac; EW0 = NW*EU_frac; IW0 = NW*IU_frac;
        yinit = [SH0; EH0; AH0; DH0; Ie0; SU0; EU0; IU0; SW0; EW0; IW0];
    end
    options = odeset('AbsTol',1e-10,'RelTol',1e-10);
    for irun = 1:50
        [~,y] = ode45(@BaseModel,linspace(0,10^4,50),yinit,options,P);
        if norm(y(end,:)-y(end-1,:))<10^-5
            break
        end
        yinit = y(end,:);
    end

    SS_mat(6,1:11) = y(end,:);

    SS_mat(6,index_m) = 1;
    SS_mat(6,index_w) = Wol_EEp_sta;
    SS_mat(6,index_all) = SS_mat(6,index_w)*SS_mat(6,index_m);    
    % 
    % if Wol_EEp_sta==1
    %     SS_mat(6,end) = 1;
    % else
    %     SS_mat(6,end) = 0;
    % end
    
    SW_frac = y(end,9)/sum(y(end,9:11),2);
    EW_frac = y(end,10)/sum(y(end,9:11),2);
    IW_frac = y(end,11)/sum(y(end,9:11),2);

end

%% (row 5) EE-EE-: malaria endemic and unstable Wolbachia endemic
NU = Wol_EEm(1); NW = Wol_EEm(2);
R0m = Cal_R0_malaria(NU,NW,P);
% if ~isnan(NU)
%     keyboard
% end
if R0m>1 && ~isnan(NU) % if Wolbachia EE- exist
    
    if flag_nearby && ~isnan(sum(SS_mat_old(5,1:11),2))
        yinit = SS_mat_old(5,1:11);
    else
        SH0 = P.gH/P.muH*0; EH0 = P.gH/P.muH*0.3; AH0 = P.gH/P.muH*0.35; DH0 = P.gH/P.muH*0.35; Ie0 = 0;
        SU0 = NU*SU_frac; EU0 = NU*EU_frac; IU0 = NU*IU_frac;
        if exist('SW_frac','var')
            SW0 = NW*SW_frac; EW0 = NW*EW_frac; IW0 = NW*IW_frac;
        else
            SW0 = NW*SU_frac; EW0 = NW*EU_frac/2; IW0 = NW*IU_frac/2;
        end
        yinit = [SH0; EH0; AH0; DH0; Ie0; SU0; EU0; IU0; SW0; EW0; IW0];
    end
    % solve a reduced system so that SU+EU+IU == NU and SW+EW+IW == NW are constant for mosquitoes
    yinit([8,11])=[]; 
    F_prop = @(x) Model_RHS(x,P,NU,NW);
    options = optimoptions('fsolve','Display','none','OptimalityTolerance', 1e-25);
    [ysol,err,~,~,~] = fsolve(F_prop,yinit,options);
    if max(max(abs(err)))>10^-5
        disp('not converged')
        keyboard
    end
    SS_mat(5,[1,2,3,4,5,6,7,9,10]) = ysol;
    SS_mat(5,8) = NU-ysol(6)-ysol(7);
    SS_mat(5,11) = NW-ysol(8)-ysol(9);

    SS_mat(5,index_m) = 1;    
    if Wol_EEm_sta~=0; keyboard;end % wolbachia EE- should always be unstable
    SS_mat(5,index_w) = Wol_EEm_sta;
    SS_mat(5,index_all) = SS_mat(5,index_w)*SS_mat(5,index_m);

    % SS_mat(5,end) = 0;

end
end

function dy = Model_RHS(y,P,NU,NW)
SH = y(1);
EH = y(2);
AH = y(3);
DH = y(4);
Ie = y(5);

SU = y(6);
EU = y(7);
IU = NU-SU-EU;

SW = y(8);
EW = y(9);
IW = NW-SW-EW;

NH = SH + EH + AH + DH;
NU = SU + EU + IU;
NW = SW + EW + IW;
NM = NU + NW;

[rho, phi, psi] = sigmoid_prob(Ie/NH, P);

BM = P.bm*P.bh*NH/(P.bm*NM+P.bh*NH);
BH = P.bm*P.bh*NM/(P.bm*NM+P.bh*NH);

LamH = BH*P.betaM*(IU+IW)/NM;
LamM = BM*(P.betaD*DH+P.betaA*AH)/NH;

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
% dIU = P.sigma*EU - P.mufu*IU;
dSW = -P.alpha*LamM*SW + gW - P.mufw*SW;
dEW = P.alpha*LamM*SW - P.sigma*EW - P.mufw*EW;
% dIW = P.sigma*EW - P.mufw*IW;

% dy = [dSH; dEH; dAH; dDH; dIe; dSU; dEU; dIU; dSW; dEW; dIW];
dy = [dSH; dEH; dAH; dDH; dIe; dSU; dEU; dSW; dEW];

end