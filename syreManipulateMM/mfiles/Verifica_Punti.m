% 07 06 07
% verifica punto di lavoro
% carica mappe
% calcola la tensione dato un vettore corrente

cd ..
load mat\fdfq_idiq_n128.mat
cd m

c_path = cd;

R = 0.022;  % 90°

% gamma ottimo @ 291A
gamma = (69:1:73) * pi/180;
Id_0 = sqrt(2) * 291 * cos(gamma);            % peak
Iq_0 = sqrt(2) * 291 * sin(gamma);            % peak
% end gamma ott

% confronto con prove berto - 08-06-07
% Id_0 = sqrt(2) * 291 * cos(gamma);            % peak
% Iq_0 = sqrt(2) * 265;            % peak


Fd =interp2(Id,Iq,Fd,Id_0,Iq_0)
Fq=interp2(Id,Iq,Fq,Id_0,Iq_0)

n = 1000; %rpm
w = 3 * pi/30 * n;
e = w * sqrt(Fd.^2 + Fq.^2) * sqrt(1/2)

I_ = Id_0 + j * Iq_0
F_ = Fd + j * Fq
E_ = j*w * F_

V_ = E_ + R * I_;

Vmod = abs(V_) * sqrt(1/2);
Imod = abs(I_) * sqrt(1/2);
gamma = angle(I_) * 180/pi;

Pelt = [];
Papp = [];
PF = [];

for jj = 1 : length(Vmod)
Pelt(jj) = 3/2*real(V_(jj) * I_(jj)');
Papp(jj) = 3/2*abs(V_(jj) * I_(jj)');
PF(jj) = Pelt(jj) / Papp(jj);
end

WP = [];

WP.kW = Pelt / 1000;
WP.Vmod = Vmod;
WP.Imod = Imod;
WP.Id = Id_0;
WP.Iq = Iq_0;
WP.gamma = gamma;

WP.PF = PF;

WP.Fd = Fd;
WP.Fq = Fq;
clc
format short
WP

% 291 A - gamma ottimo
figure
subplot(2,1,1);
pl = plot(WP.gamma,[WP.kW/250;WP.Vmod/370;WP.PF],'LineWidth',[1.5]);
grid on
legend(pl,'kW/250kW','Vph/370rms','PF')
ylim([0,2]);
title('Optimal gamma - 291 Arms','FontSize',[12],'FontWeight','Bold')
subplot(2,1,2)
pl = plot(WP.gamma,WP.Fd,'LineWidth',[1.5]);
grid on
legend(pl,'\lambda_d (Vs)');
ylim([0,2]);
xlabel('gamma (°)','FontWeight','Bold')
% 
% % prove berto 08 06 07
% K_Vs_Fref = 5146.8  % Vs to Ref Flux (control panel)
% Exp.kW = [240 239.3 241];
% Exp.Vmod = [353 359 370];
% Exp.Imod = [290 285 281];
% Exp.PF = [0.78 0.78 0.77];
% 
% Exp.Fd = [7456 7650 8000]/K_Vs_Fref;
% 
% figure
% pl = plot(WP.Fd,[WP.kW/250;WP.Vmod/370;WP.PF;WP.Imod/291],'LineWidth',[1.5]);
% grid on
% hold on
% plot(Exp.Fd,[Exp.kW/250;Exp.Vmod/370;Exp.PF;Exp.Imod/291],'x','LineWidth',[1.5]);
% legend(pl,'kW/250kW','Vph/370rms','PF','I/291rms')
% % ylim([0,2]);
% title('FEM (cont line) - ExpTest080607 (star)','FontSize',[12],'FontWeight','Bold')
% xlabel('\lambda_d (Vs)','FontWeight','Bold')
% 
% save([ c_path(1:end-1) 'mat\ExpTest_070608_250kW']);

