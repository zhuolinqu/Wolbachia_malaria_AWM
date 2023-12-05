clear all; clc
figure(1); clf
figure(2); clf
figure(3); clf

load('parameters_P_baseline.mat')
lower = [P.h_lower P.rA_lower P.rD_lower P.phiU_lower P.phiW_lower P.mufu_lower P.mufw_lower ...
    P.sigma_lower P.betaD_lower P.betaA_lower P.betaM_lower ...
    P.de_lower P.vw_lower P.ci_lower P.alpha_lower ...
    P.phis2_lower P.phir2_lower P.rhos2_lower P.rhor2_lower P.psis2_lower P.psir2_lower...
    P.cS_lower P.cE_lower P.cA_lower P.cD_lower P.dummy_lower];
upper = [P.h_upper P.rA_upper P.rD_upper P.phiU_upper P.phiW_upper P.mufu_upper P.mufw_upper ...
    P.sigma_upper P.betaD_upper P.betaA_upper P.betaM_upper ...
    P.de_upper P.vw_upper P.ci_upper P.alpha_upper ...
    P.phis2_upper P.phir2_upper P.rhos2_upper P.rhor2_upper P.psis2_upper P.psir2_upper...
    P.cS_upper P.cE_upper P.cA_upper P.cD_upper P.dummy_upper];
labs={'h','rA','rD','phiU','phiW','mufu','mufw',...
    'sigma','betaD','betaA','betaM', ...
    'de','vw','ci','alpha',...
    'phis2','phir2','rhos2','rhor2','psis2','psir2',...
    'cS','cE','cA','cD','dummy'};
load('PRCC_result_regions_100000_26.mat')

% Region 1 - 2
% Region 2 - 6
% Region 3a - 5
% Region 3b - 4
% Region 4 - 3

for k = [1 2 6 5 4 3]
    for j = 1:27
        if j <= 9
            j1 = j;
            figure(1)
        elseif j <=18
            j1 = j-9;
            figure(2)
        else
            j1 = j-18;
            figure(3)
        end
        subplot(3,3,j1)
        bin_num = 50;
        if j < 27
        if k==1
            histogram(sample_label(:,j),'binwidth',(upper(j)-lower(j))/bin_num,'binlimits',[lower(j) upper(j)],'Normalization','count')
            xlabel(labs{j})
            set(gca,'fontsize',14)
            xlim([lower(j) upper(j)])
            set(gca,'xtick',[lower(j) upper(j)])
        else
            histogram(sample_label(sample_label(:,27)==k,j),'binwidth',(upper(j)-lower(j))/bin_num,'binlimits',[lower(j) upper(j)],'Normalization','count')
        end
        else
        if k==1
            j2 = 26;
            histogram(sample_label(:,j2),'binwidth',(upper(j2)-lower(j2))/bin_num,'binlimits',[lower(j2) upper(j2)],'Normalization','count')
            xlabel(labs{j2})
            set(gca,'fontsize',14)
            xlim([lower(j2) upper(j2)])
            set(gca,'xtick',[lower(j2) upper(j2)])
        else
            histogram(sample_label(sample_label(:,27)==k,j2),'binwidth',(upper(j2)-lower(j2))/bin_num,'binlimits',[lower(j2) upper(j2)],'Normalization','count')
        end
        end
        
        hold on
    end
end
legend('All','1','2','3a','3b','4')
