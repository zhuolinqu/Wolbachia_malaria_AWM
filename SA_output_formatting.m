function [lP_list,lQ,lQ_title] = SA_output_formatting(lP_list,lQ,flag_dollar)
if flag_dollar
    lP_list = cellfun(@(x) join(["$\",x,"$"],""),lP_list,'UniformOutput',false);
end
for ip = 1:length(lP_list)
    lP_list{ip} = strrep(lP_list{ip},'\de','d_e');
    lP_list{ip} = strrep(lP_list{ip},'\rD','r_D');
    lP_list{ip} = strrep(lP_list{ip},'\rA','r_A');
    lP_list{ip} = strrep(lP_list{ip},'\mufu','\mu_{fu}');
    lP_list{ip} = strrep(lP_list{ip},'\mufw','\mu_{fw}');
    lP_list{ip} = strrep(lP_list{ip},'\phiW','\phi_w');
    lP_list{ip} = strrep(lP_list{ip},'\phiU','\phi_u');
    lP_list{ip} = strrep(lP_list{ip},'\vw','v_w');
    lP_list{ip} = strrep(lP_list{ip},'\ci','c_i');
    lP_list{ip} = strrep(lP_list{ip},'\c','c_');
    lP_list{ip} = strrep(lP_list{ip},'betaM','beta_M');
    lP_list{ip} = strrep(lP_list{ip},'betaA','beta_A');
    lP_list{ip} = strrep(lP_list{ip},'betaD','beta_D');
    lP_list{ip} = strrep(lP_list{ip},'\psir2','r_{\psi}');
    lP_list{ip} = strrep(lP_list{ip},'\psis2','s_{\psi}');
    lP_list{ip} = strrep(lP_list{ip},'\rhor2','r_{\rho}');
    lP_list{ip} = strrep(lP_list{ip},'\rhos2','s_{\rho}');
    lP_list{ip} = strrep(lP_list{ip},'\phir2','r_{\phi}');
    lP_list{ip} = strrep(lP_list{ip},'\phis2','s_{\phi}');
    lP_list{ip} = strrep(lP_list{ip},'\h','h');
    lP_list{ip} = strrep(lP_list{ip},'\dummy','dummy');
end

lQ_title = lQ;
for iq = 1:length(lQ)
    if strcmp(lQ{iq}, 'R0w');  lQ_title{iq} = '$\mathcal{R}_0^w$'; end
    if strcmp(lQ{iq}, 'R0m');  lQ_title{iq} = '$\mathcal{R}_0^m$'; end
    if strcmp(lQ{iq}, 'bifur_region');  lQ_title{iq} = 'bifur region'; end
end

end
