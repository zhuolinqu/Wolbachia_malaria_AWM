clear all; clc; figure(1); clf; figure(2); clf;
Baseline_params_malaria
P = Baseline_params_stephensi(P);
[P.ci P.vw P.alpha]

vw_vec = linspace(0.85,0.99,100);%0.85:.01:.99;
ci_vec = linspace(0.1,0.88,100);%.1:.01:.99;

outW = NaN(length(vw_vec),length(ci_vec));
outM = NaN(length(vw_vec),length(ci_vec));
outci = NaN(length(vw_vec),length(ci_vec));
outvw = NaN(length(vw_vec),length(ci_vec));
outWB = NaN(length(vw_vec),length(ci_vec));

for c1 = 1:length(vw_vec)
    P.vw = vw_vec(c1);
    P.vu = 1-P.vw;
    for c2 = 1:length(ci_vec)
        P.ci = ci_vec(c2);
        
        outci(c1,c2) = P.ci;
        outvw(c1,c2) = P.vw;
        
        % Find A point
        R0w = Cal_R0w_Del(P);
        P.phiW = R0w*P.phiU*P.mufw/(P.vw*P.mufu)+1e-8;
        SS_mat = EquilibriumState_w(P);
        
        Su = SS_mat(3,1);
        Sw = SS_mat(3,2);
        [R0m, RHM, RMH] = Cal_R0_malaria(Su,Sw,P);
        %R0m = Cal_R0m_Del(P);
        outW(c1,c2) = R0w;
        outM(c1,c2) = R0m;
        
        % Find B point, x-value
        [phiW, R0mB] = Cal_phiW_R0m(P);
        if isnan(phiW)
            %keyboard
        end
        P.phiW = phiW;
        
        [R0wB, G0w, G0u] = Cal_R0_wolbachia(P);
        outWB(c1,c2) = R0wB;
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

figure(3)
h=imagesc(ci_vec,vw_vec,outWB);
set(h,'alphadata',~isnan(outWB))
set(gca,'Ydir','normal')
title('R_0^W of B')
set(gca,'fontsize',20)
xlabel('c_i')
ylabel('v_w')
colorbar

yci = 1-(1-(1-outvw)./outvw.*outW)-outci;

figure(1)
hold on
contour(ci_vec,vw_vec,yci,[0 0],'r')
hold off

figure(2)
hold on; 
contour(ci_vec,vw_vec,outM,[1 1],'r')
hold off