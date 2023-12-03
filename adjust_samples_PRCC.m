function X = adjust_samples_PRCC(X,lP_list)
% adjust the sample when the parameter is dependent on the values from
% other paramters
[~,index] = ismember(lP_list,'phiU'); ind_phiu = find(index);
[~,index] = ismember(lP_list,'phiW'); ind_phiw = find(index);
[~,index] = ismember(lP_list,'mufu'); ind_mufu = find(index);
[~,index] = ismember(lP_list,'mufw'); ind_mufw = find(index);

[NS,~] = size(X);

for run_num=1:NS
    % convert fitness cost in lifespan -> death rate
    iP = ind_mufw; % mu_fw = mu_fu/(1-c_muf);
    X(run_num,iP) = X(run_num,ind_mufu)/(1-X(run_num,iP));
    if X(run_num,iP)<0
        keyboard
    end

    % convert fitness cost in reproduction -> egg-laying rate
    iP = ind_phiw; % phi_w = (1-c_phi)*phi_u/mu_fu * mu_fw;
    X(run_num,iP) = (1-X(run_num,iP))*X(run_num,ind_phiu)...
        /X(run_num,ind_mufu)*X(run_num,ind_mufw);
    if X(run_num,iP)<0
        keyboard
    end
end
end
