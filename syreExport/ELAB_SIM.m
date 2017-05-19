% Copyright 2015
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

% ELAB_SIM.m - Processes data from single or multiple SIMULATION from RUN_SIM

clear all, close all, clc

p = genpath('MNscripts_test_mot')
addpath(p);

load ultimo;
[FileName,pathname] = uigetfile([pathname '**.mat'],'Select the sim file');
load([pathname,FileName]);
save ultimo.mat pathname -append;

FilePath = pathname;
% valore presunto della resistenza di statore
Rs_totale=0.077;
%% %%%
MachineName=FileName(1:end-4);
temp0 = MachineName;
temp0(temp0 == '_') = '-';
temp1 = num2str(Cas.Imod,'%0.2f');
temp1(temp1 == '.') = 'A';
temp2 = num2str(Cas.gamma,'%0.1f');
temp2(temp2 == '.') = 'G';
temp3=num2str(Sim.temperature);
    
Descr= [temp0 '-' num2str(Cas.n) 'rpm-' temp1 '-' temp2 '@' temp3,'°C'];

%% Elabora i dati
%% Considera il caso di Skewing
if exist('Sim','var')
    if Sim.skew
        %% skew
        ElabSkew
        CoppiaTot = Coppia_skew1;
        TensioniTot=Va_skew;
        VaSave=Va_skew;
        FdTot = Fd_skew;
        FqTot = Fq_skew;
        
    else
        %% no skew

        ElabMov
        FdTot=Fd;
        FqTot=Fq;
        VaSave=Tensioni(:,1);
        VbSave=Tensioni(:,2);
        VcSave=Tensioni(:,3);
        
    end
else
    ElabMov
    CoppiaTot=-Coppia(:,1);
    TensioniTot=Tensioni(:,1);
    FdTot=Fd;
    FqTot=Fq;
end

% serve per dar tempo a matlab
pause(0.1)

%% FFT
n_h = round(length(Coppia360)/8);


%% Visualizza lo spettro della Coppia (360deg elt)

Ctemp = [Coppia360' Coppia(1)];
Coppia360 = repeat_n(Ctemp,round(360/Cas.sim_angle));
Coppia360 = Coppia360(1:end-1);
% 2013/05/01  MG mod: non è la posizione dell'asse d, ma l'equivalente in
% gradi del tempo di rotazione.
gradi360=theta360*180/pi-Mac.th0; % tempo preso da Magnet è in ms

%% Aromiche coppia
ArmonicheCoppia=spettro_pu(Coppia360,n_h,0);
% ArmonicheCoppia=spettro_pu(Coppia360,n_h,1000);

figure(6);
subplot(2,1,1)
plot(gradi360,Coppia360);
grid on
title(['Coppia-' Descr ' - ' datestr(now)]);

subplot(2,1,2)
bar(ArmonicheCoppia)
grid on
title(['Coppia-' Descr ' - ' datestr(now)]);
saveas(gcf,[pathname,'Coppia.fig']);

%% Visualizza lo spettro della Tensione
TensioniTot=VaSave;
Tensione = repeat_n(TensioniTot',round(360/Cas.sim_angle));

ArmonicheTensione=spettro_pu(Tensione,n_h,0);

figure(7)
subplot(2,1,1);
plot(gradi360,Tensione);
grid on
title(['Tensione-' Descr ' - ' datestr(now)]);
subplot(2,1,2);
bar(ArmonicheTensione);
grid on
title(['Tensione-' Descr ' - ' datestr(now)]);
saveas(gcf,[pathname,'Tensione.fig']);

%% Perdite nel ferro e nei magneti
if exist('pm_loss','var')
    PerditeTXT = { 'Ppm','PfeSta_h','PfeSta_c','PfeSta_o','PfeRot_h','PfeRot_c','PfeRot_o','PfeTOT','PjBarRot'};

    if ~exist('PfeSta_o','var')
        PfeSta_o = 0;
        PfeRot_o = 0;
    end
    
    if isnan(Ppm)
        Ppm=0;
    end
    PfeTOT=sum( [ Ppm , Pfes_h   , Pfes_c   , PfeSta_o , Pfer_h   , Pfer_c   , PfeRot_o, PjrBar]);
    Perdite =    [ Ppm , Pfes_h   , Pfes_c   , PfeSta_o , Pfer_h   , Pfer_c   , PfeRot_o, PfeTOT, PjrBar];
    
    figure(8);
    pie(Perdite,PerditeTXT); colormap Bone
    title(['Perdite FePM Tot= ' num2str(sum(Perdite)) ' -' Descr ' - ' datestr(now)]);
    saveas(gcf,[pathname,'PieIronLoss.fig']);

    OutDati.PM_loss=Perdite(1);
    OutDati.PfeSta_h=Perdite(2);
    OutDati.PfeSta_c=Perdite(3);
    OutDati.PfeSta_o=Perdite(4);
    OutDati.PfeRot_h=Perdite(5);
    OutDati.PfeRot_c=Perdite(6);
    OutDati.PfeRot_o=Perdite(7);
    OutDati.PjBarRot=Perdite(8);
    OutDati.Pfe=sum(Perdite(1:7));
    

end

% Joule loss
Pjs=3/2*Rs_totale*Sim.Imod(1)^2;

% 2013-05-14 modificata la valutazione della tensione concatenata e del PF,
% si tiene conto delle cadute resistive
% debug .. excel
% V_conc_pk = sqrt(3)*abs(mean(abs(FdTot+ 1i*FqTot)))*Mac.p * Cas.n*pi/30;
V_conc_pk = sqrt(3)*sqrt((Rs_totale*Id-FqTot*Mac.p * Cas.n*pi/30).^2+(Rs_totale*Iq+FdTot*Mac.p * Cas.n*pi/30).^2);
PF = cos(angle(mean(Rs_totale*Id-FqTot*Mac.p * Cas.n*pi/30) + 1i* mean(Rs_totale*Iq+FdTot*Mac.p * Cas.n*pi/30))-Cas.gamma * pi/180);
ripple_pu = std(Coppia)/mean(Coppia);

% efficiency
eta=(mean(Coppia360)*Cas.n*pi/30)/(mean(Coppia360)*Cas.n*pi/30+Ppm+Pfes_h+Pfes_c+PfeSta_o+Pfer_h+Pfer_c+PfeRot_o+PjrBar+Pjs);

Excel_out = [mean(Pjs),Sim.Imod(1),Sim.gamma(1),Sim.skew,Sim.temperature,...
    Sim.n,mean(FdTot),mean(FqTot),mean(Coppia360),ripple_pu,mean(V_conc_pk),mean(PF),mean(eta)];

out.MeanDescr={'Pjs','I' 'Gamma' 'sk' 'Temp' 'rpm' 'Fd' 'Fq' 'T' 'ripplepu' 'Vllpk' 'PF' 'eta'};
out.Mean=mean(Excel_out,1);

OutDati.Pjs=Pjs;
OutDati.Im      = mean(Excel_out(:,1),1);
OutDati.Iph     = mean(Excel_out(:,2),1);
OutDati.sk      = mean(Excel_out(:,3),1);
OutDati.Temp    = mean(Excel_out(:,4),1);
OutDati.rpm     = mean(Excel_out(:,5),1);
OutDati.Fd      = mean(Excel_out(:,6),1);
OutDati.Fq      = mean(Excel_out(:,7),1);
OutDati.T       = mean(Excel_out(:,8),1);
OutDati.ripplepu = mean(Excel_out(:,9),1);
OutDati.Vllpk   = mean(Excel_out(:,10),1);
OutDati.PF      = mean(Excel_out(:,11),1);
OutDati.eta=eta;

out.LegendaTXT=[PerditeTXT,out.MeanDescr];
out.DATA=num2cell([Perdite,out.Mean]);

StrDati=[out.LegendaTXT;out.DATA];
%% Sceda riassunto dati:
%2013/06/10 Aggiunto scheda riepilogo dati
DataLog(StrDati,MachineName,100)
saveas(gcf,[pathname,'DataLog.fig']);

save([pathname 'ElabMovHS']);


