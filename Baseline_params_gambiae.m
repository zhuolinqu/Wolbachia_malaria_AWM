function P = Baseline_params_gambiae(P)
%% Mosquito parameters taken from Childs, Hughes and Blackwood

mufu = 1/14; % death rate for adult females uninfected
mufw = 1/14; % death rate for adult females with Wolbachia
psi = 1/12; % per capita development rate
sig = 1; % per capita mating rate;
mua = 1/28; % death rate for juveniles
Ka = 100000; % carrying capacity of aquatic stage

mufup = psi/(psi+mufu)*mufu;
mufwp = psi/(psi+mufw)*mufw;

phiu = 60; % per capita egg laying rate, uninfected
phiup = psi/(psi+mua)*(sig + mufup)/(sig + mufu)*mufup/mufu*phiu;
phiupp = sig/(sig+mufup)*phiup;

phiw = 60; % per capita egg laying rate, Wolbachia infected
phiwp = psi/(psi+mua)*(sig + mufwp)/(sig + mufw)*mufwp/mufw*phiw;
phiwpp = sig/(sig+mufwp)*phiwp;

P.vw = 1; % [SA]
P.vu = 1-P.vw;
P.bf = 0.5;
P.Kf = P.bf*psi/mufup*Ka;
P.ci = 0; % new parameter

P.phiU = phiupp; % [SA]
P.phiW = phiwpp; % [SA] 
P.mufu = mufup;
P.mufw = mufwp;

%% New parameters linking malaria and Wolbachia
P.sigma = 1/10; % incubation period of malaria in mosquito
P.alpha = 0; % effectiveness of Wolbachia

end

