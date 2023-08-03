function P = Baseline_params_stephensi(P)
%% Mosquito parameters from Florez et al paper 
% psi_raw = 1/10; % from larave to adult 10 days 
% delta_raw = 1/3; % from egg to larvae 3 days 
% mu_eu_raw = 0.12;
% mu_ew_raw = 0.33;
% mu_l_raw  = 0.01;

% psi = 1/13; % from egg to larvae to adults = 1/(10+3)

% temp = delta_raw/(delta_raw+mu_eu_raw)*psi_raw/(psi_raw+mu_l_raw); % survive from egg to larave to adults  
% mu_au = psi*(1/temp -1);
% 
% temp = delta_raw/(delta_raw+mu_ew_raw)*psi_raw/(psi_raw+mu_l_raw); % survive from egg to larave to adults  
% mu_aw = psi*(1/temp -1);

P.mufu = 1/15; % Death rate for females Fu
P.mufw = 1/15; % Death rate for females Fw
 
% P.mu_mu = 1/7;  % Death rate for males Mu
% P.mu_mw = 1/7;  % Death rate for males Mw

P.phiU = 3.8; % Per capita egg-laying rate
P.phiW = 3.3; % Per capita egg-laying rate

P.vw = 1;
P.vu = 1-P.vw;
P.bf = 0.5;
P.Kf = 3e5;
P.ci = 1; % new parameter


%% New parameters linking malaria and Wolbachia
P.sigma = 1/10;
P.alpha = 0;

end