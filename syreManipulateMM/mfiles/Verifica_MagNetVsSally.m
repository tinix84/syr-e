% 12 06 07
% verifica coerenza FEM MagNet e Sally
% carica mappe
% plotta le due mappe

clear
close all

cd ..
load mat\CurveEstreme.mat
Id_MN = Id;
Iq_MN = Iq;
load mat\CurveEstreme_sally.mat
load mat\fdfq_idiq_n128.mat
cd m

clc
whos

Fd_Iq589 = interp2(Id,Iq,Fd,Id(1,:),589 * ones(1,128));

Fq_Id185 = interp2(Id,Iq,Fq,185*ones(1,128),Iq(:,1)');
% break


figure

subplot(2,1,1)
% title('SicmeRyl - FEM model')
MN = plot(Id_MN,FluxDo);
hold on
MN589 = plot(Id(1,:),Fd_Iq589,'m')
SL = plot(Id_sally,FluxD_sally,'r');
grid on
xlabel('Id - Apk')
ylabel('\lambda_d - Vs pk')
legend([MN,MN589,SL],'MagNet','MagNet,Iq = 589A','SallyMOV')

subplot(2,1,2)
MN = plot(Iq_MN,FluxQmax);
hold on
MN185 = plot(Iq(:,1),Fq_Id185,'m')
SL = plot(Iq_sally,FluxQ_sally,'r');
grid on
xlabel('Iq - Apk')
ylabel('\lambda_q - Vs pk')
legend([MN,MN185,SL],'MagNet','MagNet,Id = 185A','SallyMOV')

% MagNet integration - curve estreme.
% occhio ai nomi    - FluxQmax  è Id = 0
%                   - FluxQo    è Id = max
Fd0_model = fit(Id_MN,FluxDo,'splineinterp');
Fdmax_model = fit(Id_MN,FluxDmax,'splineinterp');

Fq0_model = fit(Iq_MN,FluxQo,'splineinterp');
Fqmax_model = fit(Iq_MN,FluxQmax,'splineinterp');

int_Fd0 = integrate(Fd0_model, Id_MN, 0);
int_Fdmax = integrate(Fdmax_model, Id_MN, 0);
int_Fq0 = integrate(Fq0_model, Iq_MN, 0);
int_Fqmax = integrate(Fqmax_model, Iq_MN, 0);

figure
plot(Id_MN,[int_Fd0 int_Fdmax])
grid on
hold on
plot(Iq_MN,[int_Fq0 int_Fqmax])
xlabel('I_{d,q} - Apk')
text(500,500,['\Delta Int - d = ' num2str(int_Fd0(end) - int_Fdmax(end))])
text(500,400,['\Delta Int - q = ' num2str(int_Fqmax(end) - int_Fq0(end))])
title('Verifica MagNet - Integrale delle Curve Estreme - Id,Iq = 0 -> 823 Apk.')