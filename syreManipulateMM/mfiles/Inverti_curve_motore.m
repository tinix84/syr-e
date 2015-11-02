%% inverti curve motore - 29 aprile 2011
% SyR style coordinates to PM style coordinates


clear all, close all, clc
addpath D:\Matlab_fdfq_idiqs

clear all, here = cd;
default_dir = 'd:\Lavoro\AOA\Database Motori';
[FILENAME, PATHNAME, FILTERINDEX] = uigetfile([default_dir '\fdfq_idiq*'], 'CARICA DATI');
run([PATHNAME 'ReadParameters'])
load([PATHNAME FILENAME])
 
% dq current range (in PM coordinates)
id_tab_min = -5;
id_tab_max = 0;
iq_tab_min = 0;
% iq_tab_min = 0;
iq_tab_max = 5;

% change the axis
tempd = Fd; tempq = Fq;
Fd = -tempq; Fq = tempd;

tempd = Id; tempq = Iq;
Id = -tempq; Iq = tempd;

% interpolazione
m = 256;    %numero righe
n = 256;   %numero colonne

[idd,iqd]=meshgrid(linspace(id_tab_min,id_tab_max,n),linspace(iq_tab_min,iq_tab_max,m));
[idq,iqq]=meshgrid(linspace(id_tab_min,id_tab_max,m),linspace(iq_tab_min,iq_tab_max,n));

fd=griddata(Id,Iq,Fd,idd,iqd);
fq=griddata(Id,Iq,Fq,idq,iqq);
% fq=fq';

%visulizzazione dati modello
figure(5),
% plot(Id(1,:),Fd'), 
% plot(Iq(:,1),Fq),
plot(idd',fd'), grid on, hold on
plot(iqq,fq),
xlabel('i_d, i_q [A]'), ylabel('\lambda_d \lambda_q [Vs]'), title(motor_name)


Id = idd;
Iq = iqq;
Fd = fd;
Fq = fq;

save([PATHNAME FILENAME(1:end-4) '_PMstyle'],'Id','Iq','Fd','Fq')
