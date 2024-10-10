function y0 = Wol_release(y0,p,yEE)
% SW EW IW % release W-infected mosquitoes
% release mosquitoes relative the yEE (w/o pre-release mosquito control)
y_release = (yEE(6)+yEE(7)+yEE(8))*p/(1-p); 
y0(9) = y_release; 
y0(10) = 0;
y0(11) = 0;
end