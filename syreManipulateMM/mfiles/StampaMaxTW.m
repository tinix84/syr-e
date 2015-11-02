
% stampa i risulati di MaxTw
[n_freq,n_coppie]=size(DatiOpt.IdMin);

pas_f=1; pas_t=1;

% torque [Nm]
T=3/2*p*(Fd.*Iq-Fq.*Id);

% Impostazione delle matrici per la visualizzazione
% delle curve iso-efficienza
[TP,VV]=meshgrid(DatiOpt.Tmap,DatiOpt.velmec);
for n=1:n_freq,
    if (DatiOpt.velmec(n)>d.velbase)
        TP(n,:)=TP(n,:)*d.velbase/DatiOpt.velmec(n); %*d.Tbase/d.TMAX;
    end
end

% stampa le isocoppie
h=figure;
[cc,hh]=contour(DatiOpt.IDtot,DatiOpt.IQtot',T,'k')%,DatiOpt.Tmap)
clabel(cc,hh); grid on
Tmax=d.TMAX;
Imax=max(max(max(Id)),max(max(Iq)));
hold on
plot((DatiOpt.IdMin(1:pas_f:end,1:pas_t:end))',(DatiOpt.IqMin(1:pas_f:end,1:pas_t:end))')
hold off
xlabel('i_d [A]'), ylabel('i_q [A]')
title(['Maximum Efficiency Locus @ Vbus= ' num2str(d.Vbus) ' V'])
axis equal, grid on
adapt_figure_fonts('Times New Roman',14,10)

% T_top_W: coppia limite
figure(2)
plot(DatiOpt.velmec(1:n_freq),DatiOpt.T_top_W,'bx-')
axis(1.1*[0 d.velmax 0 Tmax])
grid on
xlabel('Speed [rpm]')
ylabel('Torque [Nm]')
title ('Maximum Torque')

figure
[cc,hh]=contourf(VV,TP,DatiOpt.EffMOT,[75:5:85 86:1:100]);
clabel(cc,hh);
colorbar;
grid on
xlabel('Speed [rpm]')
ylabel('Torque [Nm]')
title('Motor Efficiency [%]');

figure
plot(DatiOpt.velmec(1:n_freq),DatiOpt.T_top_W,'bx-','LineWidth',1)
axis(1.1*[0 d.velmax 0 Tmax]), hold on
[cc,hh]=contour(VV,TP,DatiOpt.P_min,[0:200:1000 1500:500:4000 5000:1000:8000]);
clabel(cc,hh);
colorbar;
grid on
xlabel('Speed [rpm]')
ylabel('Torque [Nm]')
title('Motor Total Loss - W');

% figure
% % [cc,hh]=contourf(VV,TP,EffSYS,[70:1:max(max(EffSYS))]);
%
% [cc,hh]=contourf(VV,TP,DatiOpt.EffSYS,[75:5:85 86:1:95]);
% clabel(cc,hh);
% colorbar;
% grid on
% xlabel('Speed [rpm]')
% ylabel('Torque [Nm]')
% title('System Efficiency [%]');

Imin = DatiOpt.IdMin + j*DatiOpt.IqMin;
Imin = abs(Imin);

figure
plot(DatiOpt.velmec(1:n_freq),DatiOpt.T_top_W,'bx-')
axis(1.1*[0 d.velmax 0 Tmax]), hold on
%[c,h] = contour(VV,TP,Imin / sqrt(2), 50:10:max(max(Imin/sqrt(2)))); clabel(c,h)
[c,h] = contour(VV,TP,Imin , 20:20:max(max(Imin))); clabel(c,h)
grid on, hold off
xlabel('Speed [rpm]'), ylabel('Torque [Nm]')
%title('Current Map [A rms]');
title('Current Map [A pk]');

figure
plot(DatiOpt.velmec(1:n_freq),DatiOpt.T_top_W,'bx-')
axis(1.1*[0 d.velmax 0 Tmax]), hold on
[c,h] = contour(VV,TP,DatiOpt.VoMin,30:30:300); clabel(c,h)
grid on
xlabel('Speed [rpm]'), ylabel('Torque [Nm]')
title('Phase-phase motor voltage map [V pk]');

% CosfiMin
figure
plot(DatiOpt.velmec(1:n_freq),DatiOpt.T_top_W,'bx-')
axis(1.1*[0 d.velmax 0 Tmax]), hold on
[c,h] = contourf(VV,TP,DatiOpt.CosfiMin,0.7:0.05:1); clabel(c,h), colorbar
grid on
xlabel('Speed [rpm]'), ylabel('Torque [Nm]')
title('Power factor');

H_last = gcf;





