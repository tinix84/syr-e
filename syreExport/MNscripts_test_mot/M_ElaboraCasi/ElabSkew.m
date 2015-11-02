%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recupero i dati da struct Ft per ogi simulaizone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tempo espresso in [ms]
% il tempo deve essere uguale in ogni simulazione!
Tempo=Ft.values(:,1);

% Altre variabili di interesse
Idd = [];
Iqq = [];
Fdd = [];
Fqq = [];
Coppia = [];

for jj=1:Sim.skew_Nstep
      
    load([FilePath FileName(1:end-5) num2str(jj) '.mat']);
    
    ElabMov
    
    Idd(:,jj)      =   Id;
    Iqq(:,jj)      =   Iq;
    Fdd(:,jj)      =   Fd;
    Fqq(:,jj)      =   Fq;
    Va(:,jj)       =   Tensioni(:,1);
    Tcc1(:,jj)      =  -Coppia(:,1);
    Tcc2(:,jj)      =  Coppia(:,2);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elaborazione delle Misure %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tcc1 = Tcc1 - ones(size(Tcc1,1),1) * mean(Tcc1,1);
% andamento nel tempo della coppia skewata
Coppia_skew1 = mean(Tcc1,2);
Coppia_skew2 = mean(Tcc2,2);

% andamento nel tempo della Tensine Va skewata ????
Va_skew= mean(Va,2);

% Una fetta simulata con fase di corrente gamma+delta
% -> i flussi Fd Fq ottenuti vanno ruotati in avanti di delta, ovvero *exp(-i*delta)
clear i

% Matrice di angoli di rotazione
MatrAng=ones(Cas.Nstep,1)*(Cas.gamma_-Cas.gamma*pi/180)';

% Rotazione flussi su riferimento d0,q0:
temp = (Fdd + i*Fqq).*exp(-i*MatrAng);

Fd_step_skew = real(temp);
Fq_step_skew = imag(temp);
Fd_skew = mean(Fd_step_skew,2);
Fq_skew = mean(Fq_step_skew,2);

% Rotazione Correnti su riferimento d0,q0:
temp = (Idd + i*Iqq).*exp(-i*MatrAng);

Id_step_skew = real(temp);
Iq_step_skew = imag(temp);
Id_skew = mean(Id_step_skew,2);
Iq_skew = mean(Iq_step_skew,2);

CoppiaCalc_skew = 3/2*Mac.p*(Fd_step_skew.*Iq_step_skew-Fq_step_skew.*Id_step_skew);

%% debug: valori medi per ciascuna fetta
temp_Fd = mean(Fd_step_skew,1);
temp_Fq = mean(Fq_step_skew,1);
temp_Coppia = mean(Tcc1,1);
temp_V_conc_pk = sqrt(3)*abs(temp_Fd+ j*temp_Fq)*Mac.p*Cas.n*pi/30;
temp_PF = sin(Cas.gamma_' - angle(temp_Fd + j* temp_Fq))

temp_fette = [temp_Fd(3),temp_Fq(3),temp_Coppia(3),temp_V_conc_pk(3) temp_PF(3)]
temp_fette = [temp_Fd',temp_Fq',temp_Coppia',temp_V_conc_pk',temp_PF']

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Stampa Elaborazione   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on
plot(tempo, Coppia_skew1,'b','LineWidth',[2])
plot(tempo, Coppia_skew2,'r','LineWidth',[2])
plot(tempo, mean(CoppiaCalc_skew,2),'k--','LineWidth',[2])
grid on
% xlabel('tempo [s]')
% ylabel('[Nm]')
% title(['Coppia' Descr])
legend('skewed','skewed, calc','not skewed')


figure(4)
hold on
pl1 = plot(tempo,Fd_skew,'b-',tempo,Fq_skew,'r-','LineWidth',[2])
pl2 = plot(tempo,Fd_step_skew,'b--',tempo,Fq_step_skew,'r--')
grid on
% xlabel('tempo [s]')
% ylabel('[Vs]')
% title(['Flussi D(blu) e Q(rosso)' Descr])
legend([pl1(1);pl2(1)],'total','slice')



