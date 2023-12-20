%% plot wolbachia & malaria prevalence
colour_mat1 = [0 0.4470 0.7410];
colour_mat2 = [0.8500 0.3250 0.0980];
colour_mat3 = [0.9290 0.6940 0.1250];
colour_mat4 = [0.4940 0.1840 0.5560];

figure_setups; hold on
p_w = NW./(NW+NU); % wolbachia prevalence
p_m_human = (AH+DH)./sum(y(:,1:4),2); % malaria prevalence in human
p_m_mosq = (EU+IU+EW+IW)./(NW+NU);
plot(t,p_w,'-','Color',colour_mat1,'DisplayName','Wolbachia prev.')
plot(t,p_m_human,'-.','Color',colour_mat2,'DisplayName','Malaria prev. in humans') % $(A_H+D_H)$
% plot(t,p_m_human,'--','Color',colour_mat3,'DisplayName','Malaria prev. in humans') % $(A_H+D_H)$
% plot(t,p_m_human,':','Color',colour_mat4,'DisplayName','Malaria prev. in humans') % $(A_H+D_H)$
plot(t,p_m_mosq,'--','DisplayName','Malaria prev. in mosquitoes') %  $(E_M+I_M)$
xlabel('Time, days')
ylabel('Fraction of infection')
legend('Location' ,'east')
%%
% scatter(tinit+tconti_pre,0.2,300,'ko','LineWidth',4)
% scatter(222,0.33,350,'kd','filled') % peak of malaria rebound