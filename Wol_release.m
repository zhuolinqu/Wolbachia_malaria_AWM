function y0 = Wol_release(y0,p)
% SW EW IW % release W-infected mosquitoes
% disp(['release Wolbachia-carrying mosquitoes, such that there is ',num2str(p*100),' % prevalence'])
y0(9) = (y0(6)+y0(7)+y0(8))*p/(1-p); 
y0(10) = 0;
y0(11) = 0;
end