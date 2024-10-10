clearvars;
close all;
clc

%% Parameters
Baseline_params_malaria;
P = Baseline_params_stephensi(P);
P.vw = 0.95; P.vu = 1- P.vw;
tinit = 100;
tconti_pre_grid = 0:1:30;
tconti_release_grid = 1000-tinit-tconti_pre_grid;
%% Baseline scenario
tspan = 0:tinit;
SS_matM = EquilibriumState_m(P);
ySS = SS_matM(4,1:11);
options = odeset('AbsTol',1e-10,'RelTol',1e-10);
[tEE,yEE] = ode45(@BaseModel,tspan,ySS,options,P);
p_m_human_peak = NaN(size(tconti_pre_grid));
p_m_human_0 = p_m_human_peak;
for itime = 1:length(tconti_pre_grid)
    disp([num2str(itime),'/',num2str(length(tconti_pre_grid))])
    tconti_pre = tconti_pre_grid(itime);
    tconti_release = tconti_release_grid(itime);
    %% Pre-release intervention
    t = tEE;
    y = yEE;
    tspan = t(end):t(end)+tconti_pre;
    y0 = y(end,:);
    eff_malaria = 0.8; % efficacy = reduce malaria infection by X %
    y0 = malaria_control(y0,eff_malaria);
    eff_mosquito = 0;  % efficacy = reduce wild mosquito pop by X %
    y0 = mosquito_control(y0,eff_mosquito);
    if tconti_pre>0
        [tt,yy] = ode45(@BaseModel,tspan,y0,options,P);
        t = [t(1:end-1);tt];
        y = [y(1:end-1,:);yy];
        y0 = y(end,:);
    end
    %% Wolbachia release
    tspan =  t(end): t(end)+tconti_release;
    % SS_matW = EquilibriumState_w(P);
    % thres is around 49.82% = SS_matW(2,2)/(SS_matW(2,1)+SS_matW(2,2))*0.95; % Wolbachia threshold level
    p = 0.5;
    y0 = Wol_release(y0,p);
    p_m_human_0(itime) = sum(y0(3:4))/sum(y0(1:4)); % malaria prev. at Wolbachia release
    [tt,yy] = ode45(@BaseModel,tspan,y0,options,P);
    t = [t(1:end-1);tt];
    y = [y(1:end-1,:);yy];

    % Calculate the peak of malaria bounce back
    p_m_human_t = sum(yy(:,3:4),2)./sum(yy(:,1:4),2); % malaria prevalence in human (time after Wolbachia release)
    p_m_human_peak(itime) = max(p_m_human_t);
    
end
%% Plotting
figure_setups; hold on
plot(tconti_pre_grid,p_m_human_peak,'-','DisplayName','Peak of malaria following release')
plot(tconti_pre_grid,p_m_human_peak-p_m_human_0,'--','DisplayName','Rebound of malaria')
xlabel('Time since Wolbachia release, days')
ylabel('Fraction of infection')
legend('Location','best')
xlim([0 30])
% saveas(gcf,['Results/release_malaria_pre2.eps'],'epsc')

