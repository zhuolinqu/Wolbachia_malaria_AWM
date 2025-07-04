%% global SA using PRCC

% close all
clearvars
clc
format long
flag_save = 0;

% SA setting 
P.flag_adjust = 0; % if we want to adjust samples to satisfy constraints
lQ = {'R0w', 'R0m','bifur_region'};  % 'R0w', 'R0m', 'bifur_region'
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
    param = lP_list{iP};
    if P.flag_adjust == 1
        if strcmp(param,'phiW'); param = 'c_phi'; end
        if strcmp(param,'mufw'); param = 'c_muf'; end
    end
    pmin(iP,1) = P.([param,'_lower']);
    pmax(iP,1) = P.([param,'_upper']);
    pmean(iP,1) = P.(param);
end

% PRCC config
NS = 100000; % number of samples, min = k+1, 100~1000
k = length(lP_list); % # of POIs
% Pre-allocation
Size_timepts = length(time_points); % # of time points to check QOI value;
Y = NaN(NS,Size_timepts,Size_QOI);  % For each model evaluation, the QOI has dimension [Size_timepts, Size_QOI]

%% Generate parameter samples, stored in matrix X
direc = 'Results/SA/SA_PRCC/';
sample_file_name = [direc,'PRCC_sample_',num2str(NS),'_',num2str(k),'.mat'];
if ~exist(sample_file_name,'file')
    disp('generate raw samples...')
    LHS_raw = lhsdesign(NS,k); % uniform random draw with LHS sampling in [0,1]
    save(sample_file_name,'LHS_raw')
    disp('generate raw samples...DONE')
else
    disp('load raw samples...')
    load(sample_file_name,'LHS_raw')
end
label_dist = 'unif'; % 'unif', 'triangular'
X = parameterdist(LHS_raw,pmax,pmin,pmean,1,'unif'); % 'unif' 'triangular'
disp('parameter transform DONE...')
if P.flag_adjust == 1
    X = adjust_samples_PRCC(X,lP_list); % adjust samples so that it fits the biological assumptions on fitness cost
end
if flag_save; save([direc,'parameters_P_baseline.mat'],'P'); end
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
sample_label = [X,Y(:,1,strcmp('bifur_region',lQ))];
if flag_save; save([direc,'PRCC_result_Ymat_',num2str(NS),'_',num2str(k),'.mat'],'Y','P','lP_list'); end
if flag_save; save([direc,'PRCC_result_regions_',num2str(NS),'_',num2str(k),'.mat'],'Y','P','lP_list','sample_label'); end
%% PRCC on output matrix Y
% load([direc,'PRCC_result_Ymat_',num2str(NS),'_',num2str(k),'.mat'],'Y','P','lP_list','sample_label')
PRCC = NaN(k,Size_timepts,Size_QOI); stat_p = PRCC;
for itime = 1:Size_timepts
    for iQOI = 1:Size_QOI
        [rho,p] = partialcorr([X Y(:,itime,iQOI)],'type','Spearman','Rows','complete'); % ignore the NaN rows
        PRCC(:,itime,iQOI) = rho(1:end-1,end); % correlations between parameters and QOI
        stat_p(:,itime,iQOI) = p(1:end-1,end); % associated p-value
    end
end
if flag_save; save([direc,'PRCC_result_',num2str(NS),'_',num2str(k),'.mat'],'PRCC','stat_p','lP_list','lQ'); end
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
% load([direc,'PRCC_result_',num2str(NS),'_',num2str(k),'.mat'],'PRCC','stat_p','lP_list','lQ')
lP_order = {'h','rA','rD','betaA','betaD','phiU','phiW','mufu','mufw','vw','ci','alpha','sigma','betaM',...
    'de','phis2','phir2','rhos2','rhor2','psis2','psir2','cS','cE','cA','cD'};
[~,index] = ismember(lP_order,lP_list); index = index';
index(index==0)=[]; 
% [~,index] = sort(abs(PRCC(1:end-1,1,2)),'descend'); 
PRCC = PRCC([index;k],:,:); stat_p = stat_p(index,:,:,:);
lP_list = lP_list([index;k]);

%% PRCC plot 
X = categorical(lP_list);
X = reordercats(X,lP_list);
palpha = 0.01; % alpha for t-test
QOI_plot = 1:length(lQ); 
Size_QOI_plot = length(QOI_plot);
[lP_list_name,lQ,lQ_title] = SA_output_formatting(lP_list,lQ,1);
for iQOI = 1:2%Size_QOI_plot
    f = figure_setups; hold on
    set(f,'Position', [100, 55,1200, 700]);
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
    % adding dividers
    plot([5.5 5.5],[-1.2 1.2],'k-','LineWidth',1)
    plot([9.5 9.5],[-1.2 1.2],'k-','LineWidth',1)
    plot([12.5 12.5],[-1.2 1.2],'k-','LineWidth',1)
    plot([14.5 14.5],[-1.2 1.2],'k-','LineWidth',1)
    plot([21.5 21.5],[-1.2 1.2],'k-','LineWidth',1)
    if flag_save; saveas(gcf,[direc,'PRCC_result_',num2str(NS),'_',num2str(k),'_',lQ{QOI_plot(iQOI)},'.eps'],'epsc'); end
end

%%
% direc = 'Results/SA/SA_PRCC/';
% load([direc,'PRCC_result_Ymat_',num2str(NS),'_',num2str(k),'.mat'],'Y','P','lP_list','sample_label')
figure_setups;
% C = categorical(Y,[2,6,5,4,3],{'R1(2ss)','R2(6ss)','R3a(5ss)','R3b(4ss)','R4(3ss)'});
C = categorical(Y,[2,6,4,5,3],{'R1(2ss)','R2a(6ss)','R2b(4ss)','R3a(5ss)','R3b(3ss)'});
h = histogram(C);
ylim([0 10^5])
% title(['max=',num2str(max(Y)),', min=',num2str(min(Y))])

