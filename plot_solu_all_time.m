%% plot solutions in time
figure_setups;
subplot(2,2,1)
plot(t,[SH EH AH DH])
xlabel('Time, days')
ylabel('Population (human)')
legend('$S_H$','$E_H$','$A_H$','$D_H$')
ylim([0 P.gH/P.muH])

subplot(2,2,2)
plot(t,Ie./NH)
xlabel('Time, days')
ylabel('Immunity')
legend('$I_e$')

subplot(2,2,3)
plot(t,[SU EU IU])
xlabel('Time, days')
ylabel('Population (mosquito)')
legend('$S_U$','$E_U$','$I_U$')
ylim([0 P.Kf])

subplot(2,2,4)
plot(t,[SW EW IW])
xlabel('Time, days')
ylabel('Population (mosquito)')
legend('$S_W$','$E_W$','$I_W$')
ylim([0 P.Kf])