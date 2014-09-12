% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function plot_singt(out,klength,kturns,delta_sim_singt,newDir,filemot)
% single working point has been simulated

% th60 = out.SOL(:,1);
t60 = out.SOL(:,6)* klength;
t60 = [t60;t60(1)];
t = -repeat_n(t60',360/delta_sim_singt);

f_d = out.SOL(:,4)*klength*kturns;
fd60=[f_d;f_d(1)];
fd=repeat_n(fd60',6);
f_q = out.SOL(:,5)*klength*kturns;
fq60=[f_q;f_q(1)];
fq=repeat_n(fq60',6);
gamma = mean(atan2(-out.SOL(:,2),out.SOL(:,3))) * 180/pi;


delta = atan2(fq,fd) * 180/pi;

% 4 marzo 2013 - dubbio: assi PM style, gamma da asse q in senso
% antiorario
IPF = cosd(delta-gamma);

th = linspace(0,360,length(t));

%% plots
FontSize = 12;
% ymin = 1.05*round(min((t60))*100)/100;
% ymax = 1.05*round(max((t60))*100)/100;
% ytick = 1;

figure(1)
subplot(2,1,1)
pl = plot(th,abs(t)); grid on
title(['Mean Torque = ' num2str(mean(t))]);
set(pl,'LineWidth',[2]);
xlim([0 360]), %ylim([ymin ymax]),
set(gca,'FontName','TimesNewRoman','FontSize',FontSize);
ti = 0:60:360; set(gca,'XTick',ti);
%     ti = ymin:ytick:ymax; set(gca,'YTick',ti);
xl = xlabel('\theta - degrees'); set(xl,'Rotation',[0],'Fontsize',FontSize);
yl = ylabel('Nm');
set(yl,'Rotation',[90],'FontName','TimesNewRoman','Fontsize',FontSize);

figure(1);
subplot(2,1,2)
pl = plot(th,IPF);
title(['Mean IPF = ' num2str(mean(IPF))]);
grid on
set(pl,'LineWidth',2);
xl = xlabel('\theta - degrees'); set(xl,'Rotation',0,'FontName','TimesNewRoman','Fontsize',FontSize);
yl=ylabel('IPF');
set(yl,'Rotation',90,'FontName','TimesNewRoman','Fontsize',FontSize);
xlim([0 360]);
set(gca,'FontName','TimesNewRoman','FontSize',FontSize);
ti = 0:60:360;
set(gca,'XTick',ti);

saveas(gcf,[newDir filemot(1:end-4) '_T_gamma'])

ymin = round(( min(min(f_d),min(f_q)))*1000)/1000;
ymax = round(( max(max(f_d),max(f_q)))*1000)/1000;

if ymax>ymin
    ymax=1.2*ymax;
    ymin=0.8*ymin;
else
    ymax=1.2*ymin;
    ymin=0.8*ymax;
end

hdq=figure(2);
subplot(2,1,1);
hd1=plot(th,fd);
title(['Mean \lambda_d = ' num2str(mean(fd))]);
grid on
set(hd1,'LineWidth',2);
xl_hd=xlabel('\theta [Electrical degrees]');
set(xl_hd,'Rotation',0,'FontName','TimesNewRoman','Fontsize',FontSize); %,'FontWeight','Bold');
yl_hd=ylabel('\lambda_d [Wb]');
set(yl_hd,'Rotation',90,'FontName','TimesNewRoman','Fontsize',FontSize); %,'FontWeight','Bold');
xlim([0 360]);
set(gca,'FontSize',FontSize), %,'FontWeight','Bold');
ti = 0:60:360;
set(gca,'XTick',ti);

hq = subplot(2,1,2);
hq1 = plot(th,fq);
title(['Mean \lambda_q = ' num2str(mean(fq))]);
grid on
set(hq1,'LineWidth',2);
xl_hq=xlabel('\theta [Electrical degrees]');
set(xl_hq,'Rotation',0,'FontName','TimesNewRoman','Fontsize',FontSize); %'FontWeight','Bold');
yl_hq=ylabel('\lambda_q [Wb]');
set(yl_hq,'Rotation',90,'FontName','TimesNewRoman','Fontsize',FontSize); %,'FontWeight','Bold');
xlim([0 360]);
set(gca,'FontSize',FontSize); %,'FontWeight','Bold');
ti = 0:60:360;
set(gca,'XTick',ti);

saveas(hdq,[newDir filemot(1:end-4) '_plot_Flux']);
%     saveas(gcf,[newDir filemot(1:end-4) '_T_gamma'])

%% Torque Spectrum:
h = spettro_v2(abs(t),42);
NAME = filemot;
NAME(NAME=='_')='-';
title([NAME])
% saveas(gcf,[newDir filemot(1:end-4),'-torque_spectrum.fig'])
saveas(gcf,[newDir filemot(1:end-4) '-torque_spectrum.fig'])

