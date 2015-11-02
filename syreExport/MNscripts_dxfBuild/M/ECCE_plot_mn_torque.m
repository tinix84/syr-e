%% VERSIONE 20 11 2011
% plot_mn - 17 06 09

%% carica caso magnet
% disegna le stesse figure che otteniamo da matlab (T, Ts, Bt_By ..)

% clear all, close all
addpath E:\Matlab_functions
savepath = 'E:\Lavoro\Magnet\Progetti\EXPORT\';

% motor_name = Mac.MachineNameMn(1:5);
% nbs = eval(motor_name(1:2));
% nbr = eval(motor_name(4:5));

motor_name = 'CEPF37';
nbs = 12; nbr = 16;

title_string = ['macchina ' motor_name  ' - ' datestr(now)]; 

time = Ft.values(:,1);
theta_deg = 0:0.5:360; 

T = abs(Ft.values(:,12));
T0 = mean(T);

% 360 gradi
T_360 = repeat_n([T;T(1)]',6);

%% Torque
figure(1)
plot(theta_deg,T_360','r','LineWidth',[3]), grid on,
xlabel('\theta - deg'), ylabel('Nm');
axis([0 360 0 2]), title(title_string);
saveas(gcf,[savepath motor_name '_T'], 'fig')

FontSize = 10;

ymin = 0.5; ymax = 1.5;
ytick = 0.25;
%% per unit
figure(2)
pl = plot(theta_deg,T_360'/mean(T_360),'r-'), grid on,
set(pl,'LineWidth',[3]);
xlabel('\theta - deg'), ylabel('Nm');
axis([0 360 ymin ymax]), title(title_string);
saveas(gcf,[savepath motor_name '_Tpu'], 'fig')
set(gca,'FontSize',FontSize,'FontWeight','Bold');
set(gca,'PlotBoxAspectRatio',[1 0.3 1]);
ti = 0:60:360; set(gca,'XTick',ti);
ti = ymin:ytick:ymax; set(gca,'YTick',ti);
xl = xlabel('\theta - deg'), set(xl,'Rotation',[0],'Fontsize',[12],'FontWeight','Bold');
yl = ylabel('T/T_0'), set(yl,'Rotation',[0],'Fontsize',[12],'FontWeight','Bold');

ymax = 0.3;
ytick = 0.05;

%% spettro torque
h = spettro_pu(T_360(1:end-1),8 * nbs,3);
axis([0 60 0 ymax])
title([motor_name ' - spettro della coppia - ' title_string]),
colormap bone
% legend('stator teeth','stator yoke','rotor','Orientation','Vertical','Location','Best')
set(gca,'FontSize',FontSize,'FontWeight','Bold');
set(gca,'PlotBoxAspectRatio',[1 0.3 1]);
ti = 0:6:60; set(gca,'XTick',ti);
ti = 0:ytick:ymax; set(gca,'YTick',ti);
xl = xlabel('harmonic order'), set(xl,'Rotation',[0],'Fontsize',[12],'FontWeight','Bold');
yl = ylabel('T_h/T_0'), set(yl,'Rotation',[0],'Fontsize',[12],'FontWeight','Bold');
saveas(gcf,[savepath motor_name '_Tsp'], 'fig')

%%%%%
%% per unit torque
ymin = 0.5; ymax = 1.5;
ytick = 0.25;
figure(4), subplot(2,1,1)

pl = plot(theta_deg,T_360'/mean(T_360),'r-'), grid on,
set(pl,'LineWidth',[3]);
xlabel('\theta - deg'), ylabel('Nm');
axis([0 360 ymin ymax]), %title(title_string);
%saveas(gcf,[pathname filename(1:end-4) '_Tpu'], 'fig')
set(gca,'FontSize',FontSize,'FontWeight','Bold'); title(title_string)
% set(gca,'PlotBoxAspectRatio',[1 0.3 1]);
ti = 0:60:360; set(gca,'XTick',ti);
ti = ymin:ytick:ymax; set(gca,'YTick',ti);
xl = xlabel('\theta - degrees'), set(xl,'Rotation',[0],'Fontsize',[12]);
yl = ylabel('T / T_0'),
% set(yl,'Rotation',[0],'FontName','Arial','Fontsize',[12],'Position',[30 5/6*ymax]);
set(yl,'Rotation',[90],'FontName','Arial','Fontsize',[14],'FontWeight','Bold');

%% spettro torque
ymax = 0.05;
ytick = 0.025;
figure(4), subplot(2,1,2)
bar(h(1:8 * nbs),2,'r')
ylabel('p.u. della continua'), xlabel('0rdine Armonico'),grid,
% h = spettro_pu(T_360(1:end-1),8 * nbs,3);
axis([0 60 0 ymax])
% title([motor_name ' - spettro della coppia - ' title_string]),
colormap bone
% legend('stator teeth','stator yoke','rotor','Orientation','Vertical','Location','Best')
set(gca,'FontSize',FontSize,'FontWeight','Bold');
% set(gca,'PlotBoxAspectRatio',[1 0.3 1]);
ti = 0:6:60; set(gca,'XTick',ti);
ti = 0:ytick:ymax; set(gca,'YTick',ti);
xl = xlabel('h - harmonic order'), set(xl,'Rotation',[0],'Fontsize',[12],'FontWeight','Bold');
yl = ylabel('T_h / T_0'),
% set(yl,'Rotation',[0],'FontName','Arial','Fontsize',[12],'FontWeight','Bold','Position',[6 3/4*ymax]);
set(yl,'Rotation',[90],'FontName','Arial','Fontsize',[14],'FontWeight','Bold'); 
saveas(gcf,[savepath motor_name '_paper_format'], 'fig')

