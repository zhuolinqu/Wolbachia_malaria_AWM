function SS_mat = EquilibriumState_w(P)

% output [type, Equilibrium]
% row 1, DFE
% row 2, EE-
% row 3, EE+ (CIE when vw=1)

% [Fu, Fw] 

SS_mat = NaN(3,2);

[R0w, G0w, G0u] = Cal_R0_wolbachia(P);

a = P.vu/P.vw;
b = (1-P.ci)/R0w+(P.vu/P.vw-1);
c = (1-R0w)/R0w;

%% DFE
Fu = P.Kf*(1-1/G0u);
Fw = 0;
SS_mat(1,:) = [Fu, Fw];

%% EE-
if b^2-4*a*c<0 || R0w>1
    Fu = NaN; Fw = NaN;
elseif P.vw<1
    r = (-b - sqrt(b^2-4*a*c))/(2*a);
    Fu = P.Kf/(1+r*P.mufu/P.mufw)*(1-1/G0w);
    Fw = r*P.mufu/P.mufw*Fu;
elseif P.vw==1
    r = -c/b;
    Fu = P.Kf/(1+r*P.mufu/P.mufw)*(1-1/G0w);
    Fw = r*P.mufu/P.mufw*Fu;
end
SS_mat(2,:) = [Fu, Fw];

%% EE+
if b^2-4*a*c<0
    Fu = NaN; Fw = NaN;
elseif P.vw == 1 % CIE
    Fw = P.Kf*(1-1/G0w);
    Fu = 0;
elseif P.vw<1
    r = (-b + sqrt(b^2-4*a*c))/(2*a);
    Fu = P.Kf/(1+r*P.mufu/P.mufw)*(1-1/G0w);
    Fw = r*P.mufu/P.mufw*Fu;
end
SS_mat(3,:) = [Fu, Fw];

end