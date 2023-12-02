%% global SA using PRCC

close all
clearvars
clc
format long

% SA setting 
lQ = {'bifur_region'};  % 'R0w', 'R0m', 'bifur_label'
Size_QOI = length(lQ); % length of the QOI. Default = 1, unless it is an age distribution, or wants to test multiple QOIs at once
time_points = 1;
lP_list = {'h','rA','rD','phiU','phiW','mufu','mufw',...
    'sigma','betaD','betaA','betaM',...
    'de','vw','ci','alpha','phis2','phir2','rhos2','rhor2','psis2','psir2','cS','cE','cA','cD'};
lP_list{end+1} = 'dummy'; % add dummy to the POIs
Baseline_params_malaria;
P = Baseline_params_stephensi(P);
pmin = NaN(length(lP_list),1); pmax = pmin; pmean = pmin;
for iP = 1:length(lP_list)
    pmin(iP,1) = P.([lP_list{iP},'_lower']);
    pmax(iP,1) = P.([lP_list{iP},'_upper']);
    pmean(iP,1) = P.([lP_list{iP}]);
end

% PRCC config
NS = 1000; % number of samples, min = k+1, 100~1000
k = length(lP_list); % # of POIs
% Pre-allocation
Size_timepts = length(time_points); % # of time points to check QOI value;
Y = NaN(NS,Size_timepts,Size_QOI);  % For each model evaluation, the QOI has dimension [Size_timepts, Size_QOI]

%% Generate parameter samples, stored in matrix X
direc = 'Results/SA/SA_PRCC/';
if strcmp(lQ,'bifur_region')
    sample_file_name = [direc,'PRCC_sample_',num2str(NS),'_',num2str(k),'_unif.mat'];
else
    sample_file_name = [direc,'PRCC_sample_',num2str(NS),'_',num2str(k),'.mat'];
end
if ~exist(sample_file_name,'file')
    disp('generate parameter samples...')
    LHS_raw = lhsdesign(NS,k); % uniform random draw with LHS sampling in [0,1]
    if strcmp(lQ,'bifur_region')
        X = parameterdist(LHS_raw,pmax,pmin,pmean,1,'unif'); % 'unif' 'triangular'
    else
        X = parameterdist(LHS_raw,pmax,pmin,pmean,1,'triangular'); % 'unif' 'triangular'
    end
    save(sample_file_name,'X')
else
    disp('load parameter samples...')
    load([direc,'PRCC_sample_',num2str(NS),'_',num2str(k),'.mat'],'X')
end
save([direc,'parameters_P_baseline.mat'],'P')
tic
%% model evaluations
for run_num = 1:NS % Loop through each parameter sample
    if mod(run_num,NS/10)==0; disp([num2str(run_num/NS*100,3),' %']); end % display progress
    Baseline_params_malaria; % reset parameters to baseline
    P = Baseline_params_stephensi(P);
    for iP = 1:length(lP_list) % update parameter with sampled values
        P.(lP_list{iP}) = X(run_num,iP);
    end
    Q_val = QOI_value_SA(lQ,time_points,run_num,'PRCC',direc,P); % calculate QOI values
    Y(run_num,:,:) = Q_val;
end
% Y(NS,Size_timepts,Size_QOI,length(pmin),NR)
save([direc,'PRCC_result_Ymat_',num2str(NS),'_',num2str(k),'.mat'],'Y')
if strcmp(lQ,'bifur_region')
    sample_label = [X,Y];
    save([direc,'PRCC_result_regions_',num2str(NS),'_',num2str(k),'.mat'],'sample_label')
end
%% PRCC on output matrix Y
% load([direc,'PRCC_result_Ymat_',num2str(NS),'_',num2str(k),'.mat'],'Y')
PRCC = NaN(k,Size_timepts,Size_QOI); stat_p = PRCC;
for itime = 1:Size_timepts
    for iQOI = 1:Size_QOI
        [rho,p] = partialcorr([X Y(:,itime,iQOI)],'type','Spearman');
        PRCC(:,itime,iQOI) = rho(1:end-1,end); % correlations between parameters and QOI
        stat_p(:,itime,iQOI) = p(1:end-1,end); % associated p-value
    end
end
save([direc,'PRCC_result_',num2str(NS),'_',num2str(k),'.mat'],'PRCC','stat_p','lP_list','lQ')
toc

%% Plotting POIs vs. QOIs: check monotonic relationships
% for itime = 1:Size_timepts
%     figure_setups;
%     for j = 1:k
%         for iQOI = 1:Size_QOI
%             subplot(Size_QOI,k,(iQOI-1)*k+j)
%             plot(X(:,j),Y(:,itime,iQOI),'.')
%             xlabel(lP_list{j},'fontsize',14)
%             ylabel(lQ{iQOI},'fontsize',14)
%             set(gca,'fontsize',14)
%             xlim([pmin(j) pmax(j)])
%         end
%     end
%     sgtitle(['Check monotonic relationships (rows = QOIs, columns = POIs)'])
% %     sgtitle(['Check monotonic relationships at timepoint \#', num2str(itime),' (rows = QOIs, columns = parameters)'])
% end

%%
%% Sorting 
% load(['Results/vaccine_no/PRCC_result_',num2str(NS),'_',num2str(k),'.mat'],'PRCC','stat_p','lP_list','lQ')
[~,index] = sort(abs(PRCC(1:end-1,1,1)),'descend'); % sort using QOI #11, sort all the POIs except dummy
PRCC = PRCC([index;k],:,:); stat_p = stat_p(index,:,:,:);
lP_list = lP_list([index;k]);

%% PRCC plot 
X = categorical(lP_list);
X = reordercats(X,lP_list);
palpha = 0.05; % alpha for t-test
QOI_plot = 1:length(lQ); 
Size_QOI_plot = length(QOI_plot);
[lP_list_name,lQ,lQ_title] = SA_output_formatting(lP_list,lQ,1);
for iQOI = 1:Size_QOI_plot
    figure_setups; hold on
    b = bar(X,PRCC(:,1,QOI_plot(iQOI)));
    ylim([-1.1 1.1])
    xtips = b.XEndPoints;
    ytips = b.YEndPoints;
    ytips(PRCC(:,1,QOI_plot(iQOI))<0) = ytips(PRCC(:,1,QOI_plot(iQOI))<0)-0.15;
    labels = cell(1, length(lP_list_name));
    labels(stat_p(:,1,QOI_plot(iQOI))<palpha) = {'*'};
    text(xtips,ytips,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','fontsize',14)
    title(['QOI = ', lQ_title{QOI_plot(iQOI)}])
    xticklabels(lP_list_name)
    grid off
    saveas(gcf,[direc,'PRCC_result_',num2str(NS),'_',num2str(k),'_',lQ{QOI_plot(iQOI)},'.eps'],'epsc')
end


