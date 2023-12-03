function SS_mat = EquilibriumState_w(P)

% output [type, Equilibrium]
% row 1, DFE
% row 2, EE-
% row 3, EE+ (CIE when vw=1)

% last column SS_mat(:,end) indicates stability. 1 = stable, 0 = unstable

% [Fu, Fw] 

SS_mat = NaN(3,3);

[R0w, G0w, G0u] = Cal_R0_wolbachia(P);

a = P.vu/P.vw;
b = (1-P.ci)/R0w+(P.vu/P.vw-1);
c = (1-R0w)/R0w;

%% DFE
Fu = P.Kf*(1-1/G0u);
Fw = 0;
SS_mat(1,1:end-1) = [Fu, Fw];
if R0w<1
    SS_mat(1,end) = 1;
else
    SS_mat(1,end) = 0;
end

%% EE-
if b^2-4*a*c<0 || R0w>1
    Fu = NaN; Fw = NaN;
elseif P.vw==1 && R0w<1-P.ci
    Fu = NaN; Fw = NaN;
elseif P.vw<1
    r = min((-b - sqrt(b^2-4*a*c))/(2*a),(-b + sqrt(b^2-4*a*c))/(2*a));
    Fu = P.Kf/(1+r*P.mufu/P.mufw)*(1-1/G0w);
    Fw = r*P.mufu/P.mufw*Fu;
    SS_mat(2,end) = 0;
    if r<0
        Fu = NaN; Fw = NaN;
        SS_mat(2,end) = NaN;
    end
elseif P.vw==1
    r = -c/b;
    Fu = P.Kf/(1+r*P.mufu/P.mufw)*(1-1/G0w);
    Fw = r*P.mufu/P.mufw*Fu;
    SS_mat(2,end) = 0;
end
SS_mat(2,1:end-1) = [Fu, Fw];


%% EE+
if b^2-4*a*c<0
    Fu = NaN; Fw = NaN;
elseif P.vw == 1 % CIE
    Fw = P.Kf*(1-1/G0w);
    Fu = 0;
    if R0w>1-P.ci
        SS_mat(3,end) = 1;
    else
        SS_mat(3,end) = 0;
    end
elseif P.vw<1
    r = max((-b + sqrt(b^2-4*a*c))/(2*a),(-b - sqrt(b^2-4*a*c))/(2*a));
    Fu = P.Kf/(1+r*P.mufu/P.mufw)*(1-1/G0w);
    Fw = r*P.mufu/P.mufw*Fu;
    SS_mat(3,end) = 1;
    if r<0
        Fu = NaN; Fw = NaN;
        SS_mat(3,end) = NaN;
    end
end
SS_mat(3,1:end-1) = [Fu, Fw];

end