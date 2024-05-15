clear all; clc; close all
Baseline_params_malaria
P = Baseline_params_stephensi(P);
[P.ci P.vw P.alpha]

vw_vec = [linspace(0.8,0.9477,100) linspace(0.95,1,200)];
ci_vec = linspace(0,.98,200);

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

        [R0w, G0w, G0u] = Cal_R0_wolbachia(P);
        if G0w>1
            Su = SS_mat(3,1);
            Sw = SS_mat(3,2);
        else
            Su = SS_mat(1,1);
            Sw = 0;
        end
        [R0m, RHM, RMH] = Cal_R0_malaria(Su,Sw,P);
        outW(c1,c2) = R0w;
        if imag(R0m)~=0
            keyboard
        end
        if G0w < 1
            outM(c1,c2)=NaN;
        else
            outM(c1,c2) = R0m;
        end


        % Find B point, x-value
        [phiW, R0mB] = Cal_phiW_R0m(P);
        P.phiW = phiW;

        [R0wB, G0w, G0u] = Cal_R0_wolbachia(P);
        outWB(c1,c2) = R0wB;
    end
end
%% Plotting

%{
%% Heatmaps
f1 = figure_setups;
imagesc(ci_vec,vw_vec,outW)
set(gca,'Ydir','normal')
title('Minimal $\mathcal{R}_0^w$ for reduction in malaria')
xlabel('$\textit{Wolbachia}$ CI fraction, $c_i$')
ylabel('maternal transmission rate, $v_w$')
colorbar


f2 = figure_setups;
h=imagesc(ci_vec,vw_vec,outM);
set(h,'alphadata',~isnan(outM))
set(gca,'Ydir','normal')
title('$\mathcal{R}_0^m$ when \textit{Wolbachia} impacts malaria')
xlabel('$\textit{Wolbachia}$ CI fraction, $c_i$')
ylabel('maternal transmission rate, $v_w$')
colorbar

f3 = figure_setups;
h=imagesc(ci_vec,vw_vec,outWB);
set(h,'alphadata',~isnan(outWB))
set(gca,'Ydir','normal')
title('Minimal $\mathcal{R}_0^w$ for elimination of malaria')
xlabel('$\textit{Wolbachia}$ CI fraction, $c_i$')
ylabel('maternal transmission rate, $v_w$')
colorbar
set(gca,'ytick',vw_vec(1):.1:vw_vec(end))

yci = 1-(1-(1-outvw)./outvw.*outW)-outci;

figure(f1)
hold on
contour(ci_vec,vw_vec,yci,[0 0],'r','linewidth',3)
set(gca,'ytick',vw_vec(1):.1:vw_vec(end))
hold off

figure(f2)
hold on;
contour(ci_vec,vw_vec,outM,[1 1],'r','linewidth',3)
set(gca,'ytick',vw_vec(1):.1:vw_vec(end))
hold off

figure(f3)
hold on;
contour(ci_vec,vw_vec,outWB,[1 1],'r','linewidth',3)
hold off
%}
%% Contour plot
outWB1 = outWB;
outWB1(isnan(outWB)) = -.01;

f4 = figure_setups;
contour(ci_vec,vw_vec,yci,[0 0],'k','linewidth',3)
hold on
contour(ci_vec,vw_vec,outWB,[1 1],'k','linewidth',3)
[cc, hh] = contour(ci_vec,vw_vec,outWB1,[0 0],'k','linewidth',3);
xlabel('$\textit{Wolbachia}$ CI fraction, $c_i$')
ylabel('maternal transmission rate, $v_w$')
% ci values [.1 .9 .01 .2 .95 .97];
% vw values [.85  .9 .98 .98 .96 .99];
x1 = [.06 .5 0.003 .2 .6 .9];
y1 = [.875 .875 .945 .96 .98 .992];
labs = {'I','II','III','IV','V','VI'};
for j = 1:6
    text(x1(j),y1(j),labs{j})
end
hold off


