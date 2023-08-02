clear all; clc

Baseline_params_malaria
P = Baseline_params_gambiae(P);

%{
P.vw = .9;
P.vu = 1-P.vw;
P.alpha = 0.1;
P.ci = 1;
%}

phiW_list = P.phiW;
%phiW_list = linspace(1,20,300);
%phiW_list = P.phiW;

for k = 1:length(phiW_list)
    P.phiW = phiW_list(k);


    G0w = P.vw*P.bf*P.phiW/P.mufw;
    G0u = P.bf*P.phiU/P.mufu;
    R0w = P.vw*P.mufu*P.phiW/(P.mufw*P.phiU);

    a = P.vu/P.vw;
    b = (1-P.ci)/R0w+(P.vu/P.vw-1);
    c = (1-R0w)/R0w;

    rp = (-b + sqrt(b^2-4*a*c))/(2*a);
    rn = (-b - sqrt(b^2-4*a*c))/(2*a);

    for j = 1:3
        if j==1
            r = rp;
            Sustar = P.Kf/(1+r*P.mufu/P.mufw)*(1-1/G0w);
            Swstar = r*P.mufu/P.mufw*Sustar;

        elseif j ==2
            r = rn;
            Sustar = P.Kf/(1+r*P.mufu/P.mufw)*(1-1/G0w);
            Swstar = r*P.mufu/P.mufw*Sustar;

        elseif j==3
            Sustar = P.Kf*(1-1/G0u);
            Swstar = 0;
        end


        R0MH = P.bm*P.betaM*...
            (P.alpha*Swstar/(Sustar+Swstar)*P.sigma/(P.mufw+P.sigma)*1/P.mufw+...
            Sustar/(Sustar+Swstar)*P.sigma/(P.mufu+P.sigma)*1/P.mufu);

        R0HM = P.bh*P.h*(P.h+P.muH)*...
            (P.rho*P.betaD*1/(P.rD+P.muD+P.muH)+...
            (1-P.rho)*P.betaA*1/(P.rA+P.muH)+...
            P.rho*(1-P.phi)*P.rD/(P.rD+P.muD+P.muH)*P.betaA*1/(P.rA+P.muH));

        R0_M(k,j) = sqrt(R0MH*R0HM);
        R0_W(k) = R0w;
    end

end

disp(R0_M)

figure(1)
plot(R0_W,R0_M,'linewidth',2)
hold on
set(gca,'fontsize',18)
xlabel('R0_W')
ylabel('R_0^M')
legend('EE+','EE-','DFE')
grid on

