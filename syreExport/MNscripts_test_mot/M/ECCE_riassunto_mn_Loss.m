%% VERSIONE 20 11 2011
% plot_mn_loss - 15 07 09


clear all, close all
addpath C:\Matlab_functions

nbs = 24;
nbr = 8:2:20;

pfe_s = [];
pfe_r = [];

% zero_vect = nbr * 0;
for jj = 1:length(nbr)
    nbr(jj)
    if (exist(['mn_Loss\' num2str(nbs) num2str(nbr(jj)) '_Pfe_ST.mat'],'file')==2)
        load(['mn_Loss\' num2str(nbs) num2str(nbr(jj)) '_Pfe_ST.mat']);
        load(['mn_Loss\' num2str(nbs) num2str(nbr(jj)) '_Pfe_RT.mat']);
        pfe_s = [pfe_s Pfe_ST.W_tot_x_3(2)];
        pfe_r = [pfe_r Pfe_RT.W_tot_x_3(2)];
    else
        pfe_s = [pfe_s NaN];
        pfe_r = [pfe_r NaN];
    end
end

pfe_riass = [pfe_s;pfe_r]
% plot
y_max = max(pfe_s)+max(pfe_r);
y_max = 5000;
x_min = nbr(1) - 1;
x_max = nbr(end) + 1;

figure(1), subplot(3,1,1),
bar(nbr,pfe_s), grid on, ylim([0 y_max]), xlim([x_min x_max]);
ylabel('W'); % y_max = max(pfe_r);
title(['Stator loss - n_s = ' num2str(nbs)]);
figure(1), subplot(3,1,2),
bar(nbr,pfe_r), grid on, ylim([0 y_max]), xlim([x_min x_max]);
ylabel('W');
title(['Rotor Loss - n_s = ' num2str(nbs)]);
figure(1), subplot(3,1,3),
bar(nbr,pfe_s + pfe_r), grid on, ylim([0 y_max]), xlim([x_min x_max]);
xlabel('n_r'), ylabel('W');
title(['TOTAL - n_s = ' num2str(nbs)]);
saveas(gcf,['Pfe_ns' num2str(nbs) '.fig'])

save(['mn_Loss\pfe_riass_ns' num2str(nbs)],'pfe_riass')

y_max = 6000;

load mn_Loss\pfe_riass_ns12
figure(2), subplot(3,1,1),
bar(nbr,pfe_riass','stack'), grid on, ylim([0 y_max]), xlim([x_min x_max]);
colormap bone, ylabel('W'),
load mn_Loss\pfe_riass_ns18
figure(2), subplot(3,1,2),
bar(nbr,pfe_riass','stack'), grid on, ylim([0 y_max]), xlim([x_min x_max]);
ylabel('W');legend('stator loss','rotor loss','Orientation','Vertical','Location','Best')
load mn_Loss\pfe_riass_ns24
figure(2), subplot(3,1,3),
bar(nbr,pfe_riass','stack'), grid on, ylim([0 y_max]), xlim([x_min x_max]);
xlabel('n_r'), ylabel('W');




