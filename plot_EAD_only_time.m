%% plot diseased groups only
% figure_setups;
% subplot(1,2,1)
% plot(t,[EH AH DH])
% xlabel('Time, days')
% ylabel('Population (human)')
% legend('$E_H$','$A_H$','$D_H$')
% subplot(1,2,2)
% plot(t,[EH AH DH]./NH)
% xlabel('Time, days')
% ylabel('Proportion (human)')
% legend('$E_H$','$A_H$','$D_H$')
%%
figure_setups; hold on
plot(t,(SW+EW+IW)./NM,'-','DisplayName','Wolbachia prev. in mosquitoes')
plot(t,(AH+DH)./NH,'-.','DisplayName','Malaria prev. in humans')
plot(t,(IU+IW)./NM,'--','DisplayName','Malaria prev. in mosquitoes')
xlabel('Time, days')
ylabel('Fraction of infection')
legend
