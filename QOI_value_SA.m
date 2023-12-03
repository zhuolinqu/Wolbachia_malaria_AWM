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
            num_ss = sum(EquilibriumState_m_exist(P));
            if num_ss<2 || num_ss>6
                keyboard
                % 2, 3, 4, 5, 6 steady states
            end
            Q_val(:,iQ)  = num_ss;
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