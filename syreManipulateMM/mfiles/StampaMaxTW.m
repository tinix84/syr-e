
% stampa i risulati di MaxTw
[n_vel,n_coppie]=size(DatiOpt.IdMin);

pas_f=1; pas_t=1;


% Impostazione delle matrici per la visualizzazione
% delle curve iso-efficienza
[TP,VV]=meshgrid(DatiOpt.Tmap,DatiOpt.velmec);
% for n=1:n_vel
%     if (DatiOpt.velmec(n)>d.velbase)
%         TP(n,:)=TP(n,:)*d.velbase/DatiOpt.velmec(n); %*d.Tbase/d.TMAX;
%     end
% end

% stampa le isocoppie
figure
if ~isoctave()
    figSetting;
end
[cc,hh]=contour(DatiOpt.IDtot,DatiOpt.IQtot',T,'k');
clabel(cc,hh);
% [cc,hh]=contour(DatiOpt.IDtot,-DatiOpt.IQtot',-T,'k');
% clabel(cc,hh);
Tmax=d.Tmax;
Imax=max(max(max(Id)),max(max(Iq)));
plot((DatiOpt.IdMin(1:pas_f:end,1:pas_t:end))',(DatiOpt.IqMin(1:pas_f:end,1:pas_t:end))')
xlabel('$$i_d$$ [A]'), ylabel('$$i_q$$ [A]')
title(['Maximum Efficiency Locus @ Vbus= ' num2str(d.Vbus) ' V'])
axis equal

% T_top_W: coppia limite
figure
if ~isoctave()
    figSetting;
end
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_top_W,'bx-')
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_bot_W,'bx-')
grid on
xlabel('Speed [rpm]')
ylabel('Torque [Nm]')
title ('Torque Limit')

% Efficienza
figure
if ~isoctave()
    figSetting;
end
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_top_W,'-k')
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_bot_W,'-k')
[cc,hh]=contourf(VV,TP,DatiOpt.EffMOT,[64:2:86 87:1:100]);
clabel(cc,hh);
colorbar;
xlabel('Speed [rpm]')
ylabel('Torque [Nm]')
title('Motor Efficiency [\%]');

% Potenza
figure
if ~isoctave()
    figSetting;
end
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_top_W,'-k')
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_bot_W,'-k')
[cc,hh]=contourf(VV,TP,abs(DatiOpt.Potenza));
clabel(cc,hh);
colorbar;
xlabel('Speed [rpm]')
ylabel('Torque [Nm]')
title('Shaft Power [W]')

% Perdite totali
figure
if ~isoctave()
    figSetting;
end
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_top_W,'-k')
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_bot_W,'-k')
[cc,hh]=contourf(VV,TP,DatiOpt.P_min);%,[0:200:1000 1500:500:4000 5000:1000:8000]);
clabel(cc,hh);
colorbar;
xlabel('Speed [rpm]')
ylabel('Torque [Nm]')
title('Motor Total Loss - W');

% Perdite Joule Statore
figure
if ~isoctave()
    figSetting;
end
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_top_W,'-k')
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_bot_W,'-k')
[cc,hh]=contourf(VV,TP,DatiOpt.PERDITE_JOULE);
clabel(cc,hh);
colorbar;
xlabel('Speed [rpm]')
ylabel('Torque [Nm]')
title('Joule Stator Loss - W');

% figure
% if ~isoctave()
%     figSetting;
% end
% plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_top_W,'-k')
% plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_bot_W,'-k')
% [cc,hh]=contourf(VV,TP,DatiOpt.PERDITE_BARRE_ROTORE);
% clabel(cc,hh);
% colorbar;
% xlabel('Speed [rpm]')
% ylabel('Torque [Nm]')
% title('Joule Rotor Loss - W');

% Perdite Ferro
figure
if ~isoctave()
    figSetting;
end
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_top_W,'-k')
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_bot_W,'-k')
[cc,hh] = contourf(VV,TP,DatiOpt.PERDITE_FERRO);
clabel(cc,hh);
colorbar;
xlabel('Speed [rpm]')
ylabel('Torque [Nm]')
title('Core Loss - W');

% Corrente
figure
if ~isoctave()
    figSetting;
end
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_top_W,'-k')
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_bot_W,'-k')
[c,h] = contourf(VV,TP,DatiOpt.IoMin);
clabel(c,h)
xlabel('Speed [rpm]')
ylabel('Torque [Nm]')
title('Current Map [A pk]');

% Tensione concatenata
figure
if ~isoctave()
    figSetting;
end
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_top_W,'-k')
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_bot_W,'-k')
[c,h] = contourf(VV,TP,DatiOpt.VoMin);
clabel(c,h)
grid on
xlabel('Speed [rpm]')
ylabel('Torque [Nm]')
title('Phase-phase motor voltage map [V pk]');

% CosfiMin
figure
if ~isoctave()
    figSetting;
end
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_top_W,'-k')
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_bot_W,'-k')
[c,h] = contourf(VV,TP,abs(DatiOpt.CosfiMin),0.4:0.05:1);
clabel(c,h)
colorbar
xlabel('Speed [rpm]')
ylabel('Torque [Nm]')
title('Power factor');

% Perdite meccaniche
figure
if ~isoctave()
    figSetting;
end
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_top_W,'-k')
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_bot_W,'-k')
[c,h] = contourf(VV,TP,DatiOpt.PERDITE_MECC);
clabel(c,h)
colorbar
xlabel('Speed [rpm]')
ylabel('Torque [Nm]')
title('Mechanical Loss [W]');

% Perdite ferro + meccaniche
figure
if ~isoctave()
    figSetting;
end
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_top_W,'-k')
plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_bot_W,'-k')
[c,h] = contourf(VV,TP,DatiOpt.PERDITE_MECC+DatiOpt.PERDITE_FERRO);
clabel(c,h)
colorbar
xlabel('Speed [rpm]')
ylabel('Torque [Nm]')
title('Iron+Mechanical Loss [W]');

% % Slip
% figure
% if ~isoctave()
%     figSetting;
% end
% plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_top_W,'-k')
% plot(DatiOpt.velmec(1:n_vel),DatiOpt.T_bot_W,'-k')
% [c,h] = contourf(VV,TP,DatiOpt.SLIP,-1:0.1:1);
% clabel(c,h)
% colorbar
% xlabel('Speed [rpm]')
% ylabel('Torque [Nm]')
% title('Slip');

H_last = gcf;





