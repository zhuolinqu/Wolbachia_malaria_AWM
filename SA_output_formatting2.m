function x_label = SA_output_formatting2(POI)

if strcmp(POI, "$\phi_u$");  x_label = [{"$\phi_u$"}, {"reproduction for \textit{W}-free mosq."}]; end
if strcmp(POI, "$\phi_w$");  x_label = [{"$\phi_w$"}, {"reproduction for \textit{W}-carrying mosq."}]; end
if strcmp(POI, "$\mu_{fu}$");  x_label = [{"$\mu_{fu}$"}, {"death rate for \textit{W}-free mosq."}]; end
if strcmp(POI, "$\mu_{fw}$");  x_label = [{"$\mu_{fw}$"}, {"death rate for \textit{W}-carrying mosq."}]; end
if strcmp(POI, "$v_w$");  x_label = [{"$v_w$"}, {"maternal transmission rate"}]; end
if strcmp(POI, "$c_i$");  x_label = [{"$c_i$"}, {"cytoplasmic incompatibility"}]; end
if strcmp(POI, "$\alpha$");  x_label = [{"$\alpha$"}, {"relative infectivity"}]; end
end
