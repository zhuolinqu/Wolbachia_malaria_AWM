function  SS_mat_exist = EquilibriumState_m_exist(P)
% determine the existence of different ss for the full system

% malaria-Wolbachia

% (row 1) DFE-DFE: no malaria and no Wolbachia
% (row 2) DFE-EE-: no malaria and unstable Wolbachia endemic
% (row 3) DFE-EE+(or CIE): no malaria and stable Wolbachia endemic
% (row 4) EE-DFE: malaria endemic and no Wolbachia
% (row 5) EE-EE-: malaria endemic and unstable Wolbachia endemic
% (row 6) EE-EE+(or CIE): malaria endemic and stable Wolbachia endemic

%[SH, EH, AH, DH, Ie, SU, EU, IU, SW, EW, IW]

SS_mat_exist = zeros(6,1);
[R0w, G0w, G0u] = Cal_R0_wolbachia(P);
% if R0w>1
%     disp('R0w>1')
% end
SS_matW = EquilibriumState_w(P);
Wol_DFE = SS_matW(1,1:end-1);
Wol_EEm = SS_matW(2,1:end-1);
Wol_EEp = SS_matW(3,1:end-1);

a = P.vu/P.vw;
b = (1-P.ci)/R0w+(P.vu/P.vw-1);
c = (1-R0w)/R0w;
deltaW = b^2-4*a*c;

% Check
% if G0u<1 || G0w<1
%     keyboard
%     % mosquito extinction, check parameter range!
% end

%% (row 1) DFE-DFE: no malaria and no Wolbachia
SS_mat_exist(1,1) = 1; % always exists

%% (row 2) DFE-EE-: no malaria and unstable Wolbachia endemic
if R0w<1 && deltaW>0 && b<0
    SS_mat_exist(2,1) = 1;
    check_r_mm(a,b,c);
end


%% (row 3) DFE-EE+(or CIE): no malaria and high Wolbachia endemic
if P.vw<1 && deltaW>0
    if R0w<1 && b<0
        SS_mat_exist(3,1) = 1; % EE+
        check_r_pp(a,b,c);
    elseif R0w>1
        SS_mat_exist(3,1) = 1; % EE+
        check_r_pp(a,b,c);
    end
elseif P.vw==1
    SS_mat_exist(3,1) = 1; % CIE
    check_r_pp(a,b,c);
end

%% (row 4) EE-DFE: malaria endemic and no Wolbachia
SU = Wol_DFE(1); SW = Wol_DFE(2); % SW = 0;
R0m = Cal_R0_malaria(SU,SW,P);
if R0m>1
    SS_mat_exist(4,1) = 1;
end

%% (row 5) EE-EE-: malaria endemic and unstable Wolbachia endemic
NU = Wol_EEm(1); NW = Wol_EEm(2);
R0m = Cal_R0_malaria(NU,NW,P);
if R0m>1
    if R0w<1 && deltaW>0 && b<0
        SS_mat_exist(5,1) = 1;
        check_r_mm(a,b,c);
    end
end


%% (row 6) EE-EE+: malaria endemic and high Wolbachia endemic
NU = Wol_EEp(1); NW = Wol_EEp(2);
R0m = Cal_R0_malaria(NU,NW,P);
if R0m>1
    if P.vw<1 && deltaW>0
        if R0w<1 && b<0
            SS_mat_exist(6,1) = 1; % EE+
            check_r_pp(a,b,c);
        elseif R0w>1
            SS_mat_exist(6,1) = 1; % EE+
            check_r_pp(a,b,c);
        end
    elseif P.vw==1
        SS_mat_exist(6,1) = 1; % CIE
        check_r_pp(a,b,c);
    end
end
end

%%
function check_r_pp(a,b,c)
r = max((-b - sqrt(b^2-4*a*c))/(2*a),(-b + sqrt(b^2-4*a*c))/(2*a));
if r<0
    keyboard
    % should never be triggered
    % SS_mat_exist(x,1) = 0;
end
end

function check_r_mm(a,b,c)
r = min((-b - sqrt(b^2-4*a*c))/(2*a),(-b + sqrt(b^2-4*a*c))/(2*a));
if r<0
    keyboard
    % should never be triggered
    % SS_mat_exist(x,1) = 0;
end
end
