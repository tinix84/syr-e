%% 10 gennaio 2010
% calcolo le perdite nel ferro nel rotore e nello statore
%% 17 gennaio 2010
% - calcolo le perdite con formula esatta e non con la varianza
% - eliminata la fondamentale da Bxy statore

r_fe = geo.r_fe;
hf = geo.hf;
beta_f = geo.beta_f;
l = geo.l;
p = geo.p;
ns = geo.ns;
lt = geo.lt;
wt = geo.wt;

Bxy = ris_sim(:,7:end);
% punti del rotore
in = 1; fin = geo.nlay; Bx0 = Bxy(:,in:fin);
in = in + geo.nlay; fin = fin + geo.nlay; By0 = Bxy(:,in:fin);
in = in + geo.nlay; fin = fin + geo.nlay; Bx1 = Bxy(:,in:fin);
in = in + geo.nlay; fin = fin + geo.nlay; By1 = Bxy(:,in:fin);
in = in + geo.nlay; fin = fin + geo.nlay; Bx2 = Bxy(:,in:fin);
in = in + geo.nlay; fin = fin + geo.nlay; By2 = Bxy(:,in:fin);
%keyboard
% punti dello statore
in = fin+1; fin = in+ns/2-1; temp_Bx3 = Bxy(:,in:fin);
in = fin+1; fin = in+ns/2-1; temp_By3 = Bxy(:,in:fin);
%keyboard
% elimino fondamentale dallo statore
Bx3 = elimina_fond_Bxy_statore(temp_Bx3);
By3 = elimina_fond_Bxy_statore(temp_By3);

% 1/Ts: tempo di campionamento
rpm = 3000;
fs = rpm / 60 * ns * p * (nsim-1);
% ripple induzione ^2
dB0_2 = (diff(Bx0)*fs).^2 + (diff(By0)*fs).^2;
dB1_2 = (diff(Bx1)*fs).^2 + (diff(By1)*fs).^2;
dB2_2 = (diff(Bx2)*fs).^2 + (diff(By2)*fs).^2;
dB3_2 = (diff(Bx3)*fs).^2 + (diff(By3)*fs).^2;
% old
% dB0_2_old = var(Bx0) + var(By0);
% dB1_2_old = var(Bx1) + var(By1);
% dB2_2_old = var(Bx2) + var(By2);

% M530-50
Sigma = 2.7778e+006;
Spessore = 0.0005;
Keddy = Sigma * Spessore^2 /12;
%% Per ritornare alla vecchia geometria decommentare Vol_R
% volume guide rotore (m3)
%Vol_R = r_fe .* hf .* beta_f * pi/180 * l * 1e-9;
% volume denti statore (m3)
Vol_S = 2 * lt * wt * l * 1e-9; 

%% Matteo 1/11/2011 calcola volume con nuova geometria:
%
%keyboard
volume_rotore;
Vol_R=Atot*l*1e-9;

%
% Pfe rot: 3 punti per guida di flusso
Pfe_pu0 = Keddy * mean(dB0_2,1) .* Vol_R;
Pfe_pu1 = Keddy * mean(dB1_2,1) .* Vol_R;
Pfe_pu2 = Keddy * mean(dB2_2,1) .* Vol_R;

% Pfe sta: un punto per ciascun dente simulato
Pfe_pu3 = Keddy * mean(dB3_2,1) * Vol_S;

% rotore: i tre valori sono molto diversi, faccio la media
Pfe_R = mean([Pfe_pu0;Pfe_pu1;Pfe_pu2]);   
Pfe_S = Pfe_pu3;

% % debug
% figure(100), subplot(3,1,1)
% bar([Pfe_pu0;Pfe_pu1;Pfe_pu2]'), title('Pfe rotore - guida'), legend('punto0','punto1','punto2')
% figure(100), subplot(3,1,2)
% bar([Pfe_pu3]'), title('Pfe statore - singolo dente')
% figure(100), subplot(3,1,3)
% bar(2 * p * [sum(Pfe_S) sum(Pfe_R)]'), title('Confronto')
% 
% % debug
% figure(200), subplot(4,1,1)
% % plot(Bx0), hold on, plot(By0,'--'), grid on, hold off,
% plot(abs(Bx0 + 1i * By0)), grid on, 
% figure(200), subplot(4,1,2)
% % plot(Bx1), hold on, plot(By1,'--'), grid on, hold off, 
% plot(abs(Bx1 + 1i * By1)), grid on, 
% figure(200), subplot(4,1,3)
% % plot(Bx2), hold on, plot(By2,'--'), grid on, hold off, 
% plot(abs(Bx2 + 1i * By2)), grid on, 
% figure(200), subplot(4,1,4)
% % plot(Bx2), hold on, plot(By2,'--'), grid on, hold off, 
% plot(abs(Bx3 + 1i * By3)), grid on, 
