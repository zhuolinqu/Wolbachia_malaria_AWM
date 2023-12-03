function Q_val = QOI_value_SA(lQ,time_points,run_num,lmethod,direc,P)

Q_val = NaN(length(time_points),length(lQ));

for iQ = 1:length(lQ)
    switch lQ{iQ}
        case 'R0w'
            [R0w, ~, ~] = Cal_R0_wolbachia(P);
            Q_val(:,iQ)  = R0w;
        case 'R0m'
            SS_matW = EquilibriumState_w(P);
            Wol_EEp = SS_matW(3,1:end-1);
            NU = Wol_EEp(1); NW = Wol_EEp(2);
            R0m = Cal_R0_malaria(NU,NW,P);
            Q_val(:,iQ)  = R0m;
        case 'bifur_region'
            % %% compute using 1st approach
            % SS_mat = EquilibriumState_m(P); % no old mat for reference, may not be accurate
            % SS_matW = EquilibriumState_w(P);
            % % rr = region label 
            % num_ww = 3 - sum(isnan(SS_matW(:,1))); % # of Wolbachia ss
            % if num_ww == 1
            %     rr = 1;
            % elseif num_ww == 2
            %     rr = 4;
            % elseif num_ww == 3
            %     minf = (SS_mat(6,3)+SS_mat(6,4))./sum(SS_mat(6,1:4),2);
            %     if abs(minf)>10^-3
            %         rr = 2;
            %     else
            %         rr = 3;
            %     end
            % end
            %% compute using 2nd approach
            num_ss = sum(EquilibriumState_m_exist(P));
            if num_ss == 2
                rr2 = 1;
            elseif num_ss == 6
                rr2 = 2;
            elseif num_ss == 5
                rr2 = 3;
            elseif num_ss == 3
                rr2 = 4;
            else
                keyboard
            end
            % if rr2 ~= rr
            %     keyboard
            % end
            Q_val(:,iQ)  = rr2;

        otherwise
            keyboard
    end
end


end

function ind = age_range_ind(a,a_start,a_end)

[~,ind1] = min(abs(a-a_start*365)); % start from a_start years old
[~,ind2] = min(abs(a-a_end*365)); % end at a_end years old

ind = ind1:ind2;
end