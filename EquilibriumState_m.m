function  SS_mat = EquilibriumState_m(P,SS_mat_old)

% malaria-Wolbachia

% (row 1) DFE-DFE: no malaria and no Wolbachia
% (row 2) DFE-EE-: no malaria and unstable Wolbachia endemic
% (row 3) DFE-EE+: no malaria and stable Wolbachia endemic
% (row 4) EE-DFE: malaria endemic and no Wolbachia
% (row 5) EE-EE-: malaria endemic and unstable Wolbachia endemic
% (row 6) EE-EE+: malaria endemic and stable Wolbachia endemic

% last column SS_mat(:,end) indicates stability. 1 = stable, 0 = unstable

SS_mat = NaN(6,12);
[R0w, G0w, G0u] = Cal_R0_wolbachia(P);

if exist('SS_mat_old','var')
    % if the nearby steady state is provided, use these values as initial
    % conditions for steady state that are numerically captured
    flag_nearby = 1;
else
    flag_nearby = 0;
end

%[SH, EH, AH, DH, Ie, SU, EU, IU, SW, EW, IW]

%% Check
if G0u<1 || G0w<1
    keyboard
    % mosquito extinction, check parameter range!
end

SS_matW = EquilibriumState_w(P);
Wol_DFE = SS_matW(1,1:end-1); Wol_DFE_sta = SS_matW(1,end);
Wol_EEm = SS_matW(2,1:end-1); Wol_EEm_sta = SS_matW(2,end);
Wol_EEp = SS_matW(3,1:end-1); Wol_EEp_sta = SS_matW(3,end);

%% (row 1) DFE-DFE: no malaria and no Wolbachia
SU = Wol_DFE(1); SW = Wol_DFE(2); EU = 0; IU = 0; EW = 0; IW = 0;
R0m = Cal_R0_malaria(SU,0,P);
SH = P.gH/P.muH; EH = 0; AH = 0; DH = 0; Ie = 0;
SS_mat(1,1:end-1) = [SH, EH, AH, DH, Ie, SU, EU, IU, SW, EW, IW];

if R0m<1 && Wol_DFE_sta==1
    SS_mat(1,end) = 1;
else
    SS_mat(1,end) = 0;
end

%% (row 2) DFE-EE-: no malaria and unstable Wolbachia endemic
SU = Wol_EEm(1); SW = Wol_EEm(2);
if ~isnan(SU) % if Wolbachia EE- exist, always unstable
    EU = 0; IU = 0; EW = 0; IW = 0;
    %     R0m = Cal_R0_malaria(SU, SW, P);
    SH = P.gH/P.muH; EH = 0; AH = 0; DH = 0; Ie = 0;
    SS_mat(2,1:end-1) = [SH, EH, AH, DH, Ie, SU, EU, IU, SW, EW, IW];
    if Wol_EEm_sta~=0; keyboard;end % wolbachia EE- should always be unstable
    SS_mat(2,end) = 0;
end

%% (row 3) DFE-EE+: no malaria and high Wolbachia endemic
SU = Wol_EEp(1); SW = Wol_EEp(2);
if ~isnan(SU) % if Wolbachia EE+ exist
    EU = 0; IU = 0; EW = 0; IW = 0;
    R0m = Cal_R0_malaria(SU, SW, P);
    SH = P.gH/P.muH; EH = 0; AH = 0; DH = 0; Ie = 0;
    SS_mat(3,1:end-1) = [SH, EH, AH, DH, Ie, SU, EU, IU, SW, EW, IW];

    if R0m<1 && Wol_EEp_sta==1
        SS_mat(3,end) = 1;
    else
        SS_mat(3,end) = 0;
    end
end

%% (row 4) EE-DFE: malaria endemic and no Wolbachia
SU = Wol_DFE(1); SW = Wol_DFE(2); % SW = 0;
R0m = Cal_R0_malaria(SU,SW,P);
if R0m>1
    SH0 = P.gH/P.muH-1; EH0 = 1; AH0 = 0; DH0 = 0; Ie0 = 0;
    SU0 = SU; EU0 = 0; IU0 = 0; SW0 = 0; EW0 = 0; IW0 = 0;
    yinit = [SH0; EH0; AH0; DH0; Ie0; SU0; EU0; IU0; SW0; EW0; IW0];
    options = odeset('AbsTol',1e-10,'RelTol',1e-10);
    [~,y] = ode45(@BaseModel,linspace(0,1000,50),yinit,options,P);
    SS_mat(4,1:end-1) = y(end,:);

    SU_frac = y(end,6)/sum(y(end,6:8),2);
    EU_frac = y(end,7)/sum(y(end,6:8),2);
    IU_frac = y(end,8)/sum(y(end,6:8),2);

    if Wol_DFE_sta==1
        SS_mat(4,end) = 1;
    else
        SS_mat(4,end) = 0;
    end
end

%% (row 6) EE-EE+: malaria endemic and high Wolbachia endemic
NU = Wol_EEp(1); NW = Wol_EEp(2);
R0m = Cal_R0_malaria(NU,NW,P);
if abs(R0m-1)<10^-3
    keyboard
end
if R0m>1 && ~isnan(NU) % if Wolbachia EE+ exist
    if flag_nearby && ~isnan(sum(SS_mat_old(6,1:end-1),2))
        yinit = SS_mat_old(6,1:end-1);
    else
        SH0 = P.gH/P.muH*0; EH0 = P.gH/P.muH*0.3; AH0 = P.gH/P.muH*0.35; DH0 = P.gH/P.muH*0.35; Ie0 = 0;
        SU0 = NU*SU_frac; EU0 = NU*EU_frac; IU0 = NU*IU_frac;
        SW0 = NW*SU_frac; EW0 = NW*EU_frac; IW0 = NW*IU_frac;
        yinit = [SH0; EH0; AH0; DH0; Ie0; SU0; EU0; IU0; SW0; EW0; IW0];
    end
    options = odeset('AbsTol',1e-10,'RelTol',1e-10);
    for irun = 1:50
        [~,y] = ode45(@BaseModel,linspace(0,1000,50),yinit,options,P);
        if norm(y(end,:)-y(end-1,:))<10^-5
            break
        end
        yinit = y(end,:);
    end

    SS_mat(6,1:end-1) = y(end,:);

    if Wol_EEp_sta==1
        SS_mat(6,end) = 1;
    else
        SS_mat(6,end) = 0;
    end

    SW_frac = y(end,9)/sum(y(end,9:11),2);
    EW_frac = y(end,10)/sum(y(end,9:11),2);
    IW_frac = y(end,11)/sum(y(end,9:11),2);

end

%% (row 5) EE-EE-: malaria endemic and unstable Wolbachia endemic
NU = Wol_EEm(1); NW = Wol_EEm(2);
R0m = Cal_R0_malaria(NU,NW,P);

if R0m>1 && ~isnan(NU) % if Wolbachia EE- exist
    if flag_nearby && ~isnan(sum(SS_mat_old(5,1:end-1),2))
        yinit = SS_mat_old(5,1:end-1);
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
    options = optimoptions('fsolve','Display','none','OptimalityTolerance', 1e-25);
    F_prop = @(x) BaseModel(0,x,P);
    [ysol,err,~,~,~] = fsolve(F_prop,yinit,options);
    if max(max(abs(err)))>10^-5
        disp('not converged')
        keyboard
    end
    SS_mat(5,1:end-1) = ysol;
    if Wol_EEm_sta~=0; keyboard;end % wolbachia EE- should always be unstable
    SS_mat(5,end) = 0;

end
end
