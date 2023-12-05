Baseline_params_malaria
P = Baseline_params_stephensi(P);

figure(1)
clf

col = lines(4);

vw_vec = [0.4 0.5 0.75 0.95];
for j1 = 1:length(vw_vec)
    phiW_vec = 0:.05:10;
    P.vw = vw_vec(j1);
    P.vu = 1-P.vw;
    for j = 1:length(phiW_vec)
        P.phiW = phiW_vec(j);
        SS_mat = EquilibriumState_w(P);
        [R0w, G0w, G0u] = Cal_R0_wolbachia(P);
        figure(1)
        xlim([0 2])
        ylim([0 1])
        plot(R0w,SS_mat(:,2)./(SS_mat(:,1)+SS_mat(:,2)),'.','color',col(j1,:))
        hold on
        
    end
end
hold off
set(gca,'fontsize',20)
xlabel('R_0^w(\phi_w)')
ylabel('Fraction infectious')
title('Changes in v_w (ci = 1)')

%%
P = Baseline_params_stephensi(P);
figure(2)
clf

col = lines(4);

ci_vec = [0.05 0.5 0.75 0.95];
for j1 = 1:length(ci_vec)
    phiW_vec = 0:.05:10;
    P.ci = ci_vec(j1);
    for j = 1:length(phiW_vec)
        P.phiW = phiW_vec(j);
        SS_mat = EquilibriumState_w(P);
        [R0w, G0w, G0u] = Cal_R0_wolbachia(P);
        figure(2)
        xlim([0 2])
        ylim([0 1])
        plot(R0w,SS_mat(:,2)./(SS_mat(:,1)+SS_mat(:,2)),'.','color',col(j1,:))
        hold on
        
    end
end
hold off
set(gca,'fontsize',14)
xlabel('R_0^w(\phi_w)')
ylabel('Fraction infectious')
title('Changes in c_i (v_w = 0.95)')