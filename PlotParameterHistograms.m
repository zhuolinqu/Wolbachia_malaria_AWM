clear all; clc
close all

file_dir = 'Results/SA/SA_PRCC/';
load([file_dir,'PRCC_result_regions_100000_26'],'P','lP_list','sample_label')
lower = NaN(length(lP_list),1); upper = lower; pmean = lower;
for iP = 1:length(lP_list)
    param = lP_list{iP};
    lower(iP,1) = P.([param,'_lower']);
    upper(iP,1) = P.([param,'_upper']);
end
[labs,~,~] = SA_output_formatting(lP_list,'',1);

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
                %xlabel(labs{j2})
                set(gca,'fontsize',14)
                xlim([0 .1])
                %xlim([lower(j2) upper(j2)])
                %set(gca,'xtick',[lower(j2) upper(j2)])
                set(gca,'xtick',[],'ytick',[])
            else
                histogram(sample_label(sample_label(:,27)==k,j2),'binwidth',(upper(j2)-lower(j2))/bin_num,'binlimits',[lower(j2) upper(j2)],'Normalization','count')
            end
        end

        hold on
    end
end
legend('All','1','2','3a','3b','4')

%% Figures for manuscript
f3 = figure_setups;

j1 = 1;
for k = [ 2 6 5 4 3]
    for j = [4:7 13:14]
        subplot(2,3,j1)
        bin_num = 50;
        histogram(sample_label(sample_label(:,27)==k,j),'binwidth',(upper(j)-lower(j))/bin_num,'binlimits',[lower(j) upper(j)],'Normalization','count')
        if k==2
            xlabel(labs{j})
            xlim([lower(j) upper(j)])
            set(gca,'xtick',[lower(j) upper(j)])
            set(gca,'ytick',[0 1000 2000],'yticklabel',[0 .5 1])
        end
        hold on

        j1 = j1 + 1;
    end
    j1  = 1;
end

% legend
f4 = figure_setups;
set(f4,'Position', [100, 55, 120, 250])
for k = [ 2 6 5 4 3]
    j = 4;
    bin_num = 50;
    histogram(sample_label(sample_label(:,27)==k,j),'binwidth',(upper(j)-lower(j))/bin_num,'binlimits',[lower(j) upper(j)],'Normalization','count')
    xlim([0 .1])
    set(gca,'xtick',[],'ytick',[])
    hold on
end

lgd = legend('1','2','3a','3b','4');
title(lgd,'Region')
set(lgd,'location','east','box','off')

%%
print(f3,'-vector', '-depsc', 'SA_par_histogram.eps')
print(f4,'-vector', '-depsc', 'SA_par_histogram_legend.eps')