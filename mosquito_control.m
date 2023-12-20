function y0 = mosquito_control(y0,eff_mosquito)
% SU EU IU % mosquito intervention
% eff_mosquito = efficacy in reducing wild mosquito pop by X %
% disp(['remove mosquito population by ',num2str(eff_mosquito*100),' %...'])
y0(6) = y0(6)+eff_mosquito*y0(7)+eff_mosquito*y0(8);
y0(7) = (1-eff_mosquito)*y0(7);
y0(8) = (1-eff_mosquito)*y0(8);