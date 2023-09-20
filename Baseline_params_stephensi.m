function P = Baseline_params_stephensi(P)
%% Mosquito parameters from Florez et al paper 
mu_fu = 1/12; % Death rate for females Fu
mu_fw = 1/14; % Death rate for females Fw
phi_u = 50/12; % Per capita egg-laying rate
phi_w = 50/14; % Per capita egg-laying rate
psi_pr = 1/13; % from egg to larvae to adults = 1/(10+3) % psi = 1/10; % FIXED; was 1/18 in Florez et al paper 

psi_raw = 1/10; % from larave to adult 10 days 
delta_raw = 1/3; % from egg to larvae 3 days 
mu_eu_raw = 0.12;
mu_ew_raw = 0.33;
mu_l_raw  = psi_raw/0.8 - psi_raw;
temp = delta_raw/(delta_raw+mu_eu_raw)*psi_raw/(psi_raw+mu_l_raw); % survive from egg to larave to adults  
mu_au = psi_pr*(1/temp -1);
temp = delta_raw/(delta_raw+mu_ew_raw)*psi_raw/(psi_raw+mu_l_raw); % survive from egg to larave to adults  
mu_aw = psi_pr*(1/temp -1);

% remove aquatic stage
% mu_fu_pr = 1/(1/psi_pr+1/mu_fu);
% mu_fw_pr = 1/(1/psi_pr+1/mu_fw);
phi_u_pr = phi_u*psi_pr/(psi_pr+mu_au);
phi_w_pr = phi_w*psi_pr/(psi_pr+mu_aw);

P.phiU = phi_u_pr;
P.phiW = phi_w_pr;

P.mufu = mu_fu; % 1/ time spent in adult females, used for malaria model
P.mufw = mu_fw;

P.vw = 1;
P.vu = 1-P.vw;
P.bf = 0.5;
P.Kf = 3e5;
P.ci = 1; % new parameter

%% New parameters linking malaria and Wolbachia
P.sigma = 1/10;
P.alpha = 0;

end