function y0 = malaria_control(y0,eff_malaria)
% SH EH AH DH % malaria intervention
% eff_malaria = efficacy in reducing malaria infection by X %
% disp(['reduce malaria infection by ',num2str(eff_malaria*100),' %...'])
y0(1) = y0(1)+eff_malaria*y0(3)+eff_malaria*y0(4);
y0(2) = y0(2);
y0(3) = (1-eff_malaria)*y0(3);
y0(4) = (1-eff_malaria)*y0(4);

end