% spettro_torque - 28 05 09

%% spettro
if exist('T','var')
    T_360 = repeat_n(T,6); mean(T_360)
    h = spettro_pu(T_360,8 * nbs,2);
    title([motor_name ' - spettro della coppia - ' title_string]),
end

if exist('T1','var')
    T1_360 = repeat_n(T1,6); mean(T1_360)
    h = spettro_pu(T1_360,8 * nbs,3);
    title([motor_name ' - spettro della coppia COENERGIA - ' title_string]),
end

if exist('T2','var')
    T2_360 = repeat_n(T2,6); mean(T2_360)
    h = spettro_pu(T2_360,8 * nbs,4);
    title([motor_name ' - spettro della coppia B x dFr - ' title_string]),
end

if exist('T3','var')
    T3_360 = repeat_n(T3,6); mean(T3_360)
    h = spettro_pu(T3_360,8 * nbs,5);
    title([motor_name ' - spettro della coppia COENERGIA - ' title_string]),
end

if exist('T4','var')
    T4_360 = repeat_n(T4,6); mean(T4_360)
    h = spettro_pu(T4_360,8 * nbs,6);
    title([motor_name ' - spettro della coppia B x dFs - ' title_string]),
end