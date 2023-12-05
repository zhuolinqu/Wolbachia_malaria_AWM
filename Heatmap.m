clear all; clc; figure(1); clf; figure(2); clf;
Baseline_params_malaria
P = Baseline_params_stephensi(P);
[P.ci P.vw P.alpha]

vw_vec = 0.1:.02:1;
ci_vec = 0.1:.01:1;

outW = NaN(length(vw_vec),length(ci_vec));
outM = NaN(length(vw_vec),length(ci_vec));
outci = NaN(length(vw_vec),length(ci_vec));
outvw = NaN(length(vw_vec),length(ci_vec));

for c1 = 1:length(vw_vec)
    P.vw = vw_vec(c1);
    P.vu = 1-P.vw;
    for c2 = 1:length(ci_vec)
        P.ci = ci_vec(c2);
        R0w = Cal_R0w_Del(P);
        SS_mat = EquilibriumState_w(P);
        outtmp(c1,c2) = SS_mat(3,1);
        outci(c1,c2) = P.ci;
        outvw(c1,c2) = P.vw;
        Su = SS_mat(3,1);
        Sw = SS_mat(3,2);
        [R0m, RHM, RMH] = Cal_R0_malaria(Su,Sw,P);
        outW(c1,c2) = R0w;
        outM(c1,c2) = R0m;
    end
end
%% Plotting 
figure(1)
imagesc(ci_vec,vw_vec,outW)
set(gca,'Ydir','normal')
title('R_0^W of A')
set(gca,'fontsize',20)
xlabel('c_i')
ylabel('v_w')
%ylabel('c_i')
%xlabel('v_w')
colorbar
hold on
contour(ci_vec,vw_vec,outW,[1 1],'k')
%contour(vw_vec,ci_vec,outW',[1 1],'k')
hold off

figure(2)
h=imagesc(ci_vec,vw_vec,outM);
set(h,'alphadata',~isnan(outM))
set(gca,'Ydir','normal')
title('R_0^M of A')
set(gca,'fontsize',20)
xlabel('c_i')
ylabel('v_w')
colorbar

yci = 1-(1-(1-outvw)./outvw.*outW)-outci;

figure(1)
hold on
contour(ci_vec,vw_vec,yci,[0 0],'r')
hold off