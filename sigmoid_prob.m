%% sigmoidal function for converting immunity to probabilties
function [rho, phi, psi] = sigmoid_prob(x, P)

sigmoid_inc = @(f_0, f_1, x, s_2, r_2) f_0 + (f_1-f_0)./(1 + exp(-(x-s_2)/r_2));

rho = sigmoid_inc(P.rhof0, P.rhof1, x, P.rhos2, P.rhor2);
phi = sigmoid_inc(P.phif0, P.phif1, x, P.phis2, P.phir2);
psi = sigmoid_inc(P.psif0, P.psif1, x, P.psis2, P.psir2);

end
