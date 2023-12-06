%% plot wolbachia & malaria prevalence
figure_setups; hold on
p_w = NW./(NW+NU); % wolbachia prevalence
p_m_human = (AH+DH)./sum(y(:,1:4),2); % malaria prevalence in human
p_m_mosq = (IU+IW)./(NW+NU);
plot(t,p_w,'-','DisplayName','Prevalence of Wolbachia')
plot(t,p_m_human,'-.','DisplayName','Prevalence of malaria (human)')
plot(t,p_m_mosq,'--','DisplayName','Prevalence of malaria (mosquito)')
xlabel('Time, days')
ylabel('Fraction of infection')
legend('Location' ,'best')