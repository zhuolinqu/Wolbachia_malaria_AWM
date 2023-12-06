%% plot sigmoids
figure_setups;
[rho, phi, psi] = sigmoid_prob(Ie./NH, P);
plot(t,rho,t,phi,t,psi)
legend('rho','phi','psi')