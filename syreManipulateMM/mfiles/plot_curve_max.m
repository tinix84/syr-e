% verify_curve_Max - 12 01 08

%% tensioni e flussi da evidenziare in AOA
% default
correnti = Imax * [0.5 0.71 1];
flussi = linspace(lm,F_A,4);
coppie = Tmax * [0.1:0.1:1];

close all
% AOA - current frame
figure(1), hold on,
plot(id_KtMax,iq_KtMax,'k','LineWidth',2.5);
plot(id_KvMax,iq_KvMax,'k','LineWidth',1.5);
[c,h]=contour(Id,Iq,I_,correnti,'--k');
[c,h]=contour(Id,Iq,F_,flussi);
[c,h]=contour(Id,Iq,T_,coppie);
h = clabel(c,h,'LabelSpacing',[144]);% colorbar;
grid on, xlabel('i_d (A)'),ylabel('i_q (A)');
title('AOA - current','FontSize',[14],'FontWeight','bold');
axis equal
% xlim([0 Imax]); ylim([0 Imax]);
% xlim([0 1000]); ylim([0 Imax]);
h=gcf(); %AS
if isoctave()
    fig_name=strcat(pathname1, 'AOA-i');
    hgsave(h,[fig_name]);
else
    saveas(gcf,[pathname1 'AOA-i'])
end

% AOA - flux frame
fig1 = figure; hold on;
[c,h]=contour(FD,FQ,I,correnti,'k'); grid on;
[c,h]=contour(FD,FQ,F,flussi);
% h = clabel(c,h,'LabelSpacing',[600]); %colorbar;
[c,h]=contour(FD,FQ,T,coppie);
h = clabel(c,h,'LabelSpacing',288);% colorbar;
plot(fd_KvMax,fq_KvMax,'k','LineWidth',1.5);
plot(fd_KtMax,fq_KtMax,'k','LineWidth',1.5);
grid on, xlabel('Fd [Vs]'),ylabel('Fq [Vs]'),zlabel('b locus');
title('AOA -Flux','FontSize',[14],'FontWeight','bold');
% xlim([0 2]);ylim([-0.1 0.1]);
axis equal
set(gca,'PlotBoxAspectRatio',[2 1 1])
% set(gca,'YTick',[-0.1 -0.05 0 0.05])
% legend('Nm/A','Nm/Vs','PF/Nm')
h=gcf(); %AS
if isoctave()
    fig_name=strcat(pathname1, 'AOA-f');
    hgsave(h,[fig_name]);
else
    saveas(gcf,[pathname1 'AOA-f'])
end

% print -dpsc2 AOA_Flux

