clear all; clc
close all

flag_save = 1;
file_dir = 'Results/SA/';
load([file_dir,'SA_PRCC/PRCC_result_regions_100000_26.mat'],'P','lP_list','sample_label')
lower = NaN(length(lP_list),1); upper = lower; pmean = lower;
for iP = 1:length(lP_list)
    param = lP_list{iP};
    lower(iP,1) = P.([param,'_lower']);
    upper(iP,1) = P.([param,'_upper']);
end
[labs,~,~] = SA_output_formatting(lP_list,'',1);

% lower = [P.h_lower P.rA_lower P.rD_lower P.phiU_lower P.phiW_lower P.mufu_lower P.mufw_lower ...
%     P.sigma_lower P.betaD_lower P.betaA_lower P.betaM_lower ...
%     P.de_lower P.vw_lower P.ci_lower P.alpha_lower ...
%     P.phis2_lower P.phir2_lower P.rhos2_lower P.rhor2_lower P.psis2_lower P.psir2_lower...
%     P.cS_lower P.cE_lower P.cA_lower P.cD_lower P.dummy_lower];
% upper = [P.h_upper P.rA_upper P.rD_upper P.phiU_upper P.phiW_upper P.mufu_upper P.mufw_upper ...
%     P.sigma_upper P.betaD_upper P.betaA_upper P.betaM_upper ...
%     P.de_upper P.vw_upper P.ci_upper P.alpha_upper ...
%     P.phis2_upper P.phir2_upper P.rhos2_upper P.rhor2_upper P.psis2_upper P.psir2_upper...
%     P.cS_upper P.cE_upper P.cA_upper P.cD_upper P.dummy_upper];
% labs={'h','rA','rD','phiU','phiW','mufu','mufw',...
%     'sigma','betaD','betaA','betaM', ...
%     'de','vw','ci','alpha',...
%     'phis2','phir2','rhos2','rhor2','psis2','psir2',...
%     'cS','cE','cA','cD','dummy'};

% sample_label (1:26) sampling
% sample_label (27) labels

% Region 1 - 2
% Region 2 - 6
% Region 3a - 5
% Region 3b - 4
% Region 4 - 3

% Region 1 - 2 SS
% Region 2 (2a & 2b) - 6 & 4 SS
% Region 3 (3a & 3b) - 5 & 3 SS

for igroup = 1:3
    if igroup ==1
        k = 2;
    elseif igroup ==2
        k = [6,4];
    elseif igroup ==3
        k = [5,3];
    end
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
            if igroup==1
                histogram(sample_label(:,j),'binwidth',(upper(j)-lower(j))/bin_num,'binlimits',[lower(j) upper(j)],'Normalization','count')
                xlabel(labs{j})
                set(gca,'fontsize',14)
                xlim([lower(j) upper(j)])
                set(gca,'xtick',[lower(j) upper(j)])
            else
                ind = [find(sample_label(:,27)==k(1)); find(sample_label(:,27)==k(2))];
                histogram(sample_label(ind,j),'binwidth',(upper(j)-lower(j))/bin_num,'binlimits',[lower(j) upper(j)],'Normalization','count')
                % ylim([0 0.1])
            end
        else
            if igroup==1
                j2 = 26;
                histogram(sample_label(:,j2),'binwidth',(upper(j2)-lower(j2))/bin_num,'binlimits',[lower(j2) upper(j2)],'Normalization','count')
                xlabel(labs{j2})
                set(gca,'fontsize',14)
                xlim([0 0.1])
                xlim([lower(j2) upper(j2)])
                set(gca,'xtick',[lower(j2) upper(j2)])
            else
                ind = [find(sample_label(:,27)==k(1)); find(sample_label(:,27)==k(2))];
                histogram(sample_label(ind,j2),'binwidth',(upper(j2)-lower(j2))/bin_num,'binlimits',[lower(j2) upper(j2)],'Normalization','count')
                % ylim([0 0.1])
            end
        end

        hold on
    end
end
legend('1','2','3')

%% Figures for manuscript
labs_plot = {'phiU','phiW','mufu','mufw','vw','ci','alpha'};
for iPOI = 1:length(labs_plot)
    figure_setups;
    POI_ind = find(strcmp(labs_plot{iPOI},lP_list));
    for igroup = 1:3
        if igroup ==1
            k = 2;
        elseif igroup ==2
            k = [6,4];
        elseif igroup ==3
            k = [5,3];
        end
        j = POI_ind;
        bin_num = 50;
        ind = [];
        for ik = 1:length(k)
            ind = [ind; find(sample_label(:,27)==k(ik))];
        end
        histogram(sample_label(ind,j),'binwidth',(upper(j)-lower(j))/bin_num,'binlimits',[lower(j) upper(j)],'Normalization','count')
        xlabel(labs{j})
        xlim([lower(j) upper(j)])
        ylim([0 2000])
        yticks([0 1000 2000])
        yticklabels([0 0.01 0.02])
        % ylim([0 0.1])
        set(gca,'xtick',[round(lower(j),2) round(upper(j),2)])
        hold on
    end
    if flag_save; saveas(gcf,[file_dir,'Histogram_result_',labs_plot{iPOI},'.eps'],'epsc'); end
end

%% legend
f4 = figure_setups;
set(f4,'Position', [100, 55, 250, 350])
for k = [3 4 5]
    j = 4;
    bin_num = 50;
    histogram(sample_label(sample_label(:,27)==k,j),'binwidth',(upper(j)-lower(j))/bin_num,'binlimits',[lower(j) upper(j)],'Normalization','count')
    xlim([0 0.1])
    set(gca,'xtick',[],'ytick',[])
    hold on
end
lgd = legend('1','2','3');
title(lgd,'Region')
set(lgd,'location','east','box','off')
if flag_save; saveas(gcf,[file_dir,'Histogram_result_legend.eps'],'epsc'); end