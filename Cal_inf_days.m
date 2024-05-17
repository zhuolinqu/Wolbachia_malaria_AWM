function days = Cal_inf_days(t,AH,DH)

days = trapz(t,AH+DH);

end