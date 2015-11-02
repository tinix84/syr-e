
tempo   =   Ft.values(:,1);         %ms
Coppia = mean([-Ft.values(:,12) Ft.values(:,13)],2);
Flussi =   Ft.values(:,2:4);         %Vs
Tensioni    =   Ft.values(:,5:7);  %Vpk fase
Id  =   Ft.values(:,8);
Iq  =   Ft.values(:,9);
Fd  =   Ft.values(:,10);
Fq  =   Ft.values(:,11);
CoppiaCalc=3/2*Mac.p*(Fd.*Iq-Fq.*Id);

%% riporta a 360 deg elettrici
if round(360/Cas.sim_angle)>1
    
    %% coppia
    Ctemp = [Coppia' Coppia(1)];
    Coppia360 = repeat_n(Ctemp,round(360/Cas.sim_angle));
    Coppia360 = Coppia360(1:end-1);
    
    dtime = diff(tempo);
    dtime = dtime(1);
    time = linspace(tempo(1),round(360/Cas.sim_angle)*tempo(end),length(Coppia360));
    
    %% flussi
    theta = (Mac.th0 * pi/180 + tempo/1000 * Cas.n * pi/30 * Mac.p);
    
    temp = [Fd' Fd(1)];
    temp360 = repeat_n(temp,round(360/Cas.sim_angle));
    Fd360 = temp360(1:end-1);
    
    temp = [Fq' Fq(1)];
    temp360 = repeat_n(temp,round(360/Cas.sim_angle));
    Fq360 = temp360(1:end-1);
    
    %% tensioni
    Vdq = abc2dq(Tensioni(:,1)',Tensioni(:,2)',Tensioni(:,3)',-theta')
    Vd = Vdq(1,:)';Vq = Vdq(2,:)';
    
%     figure(100)
%     plot(theta,[Vd Vq],'--'), hold on
%     plot(theta,Tensioni)

    temp = [Vd' Vd(1)];
    temp360 = repeat_n(temp,round(360/Cas.sim_angle));
    Vd360 = temp360(1:end-1);
    
    temp = [Vq' Vq(1)];
    temp360 = repeat_n(temp,round(360/Cas.sim_angle));
    Vq360 = temp360(1:end-1);
    
    %% dq to abc
    theta360 = (Mac.th0 * pi/180 + time/1000 * Cas.n * pi/30 * Mac.p);
    
    Flussi360 = dq2abc(Fd360,Fq360,theta360);
    Tensioni360 = dq2abc(Vd360,Vq360,theta360);
    
    %% plotto le nuove curve
    tempo = time';
    Coppia = Coppia360';
    Flussi = Flussi360';
    Fd = Fd360';
    Fq = Fq360';
    Tensioni = Tensioni360';
else
    
    Coppia360 = Coppia;
    Flussi360 = Flussi;
    Tensioni360 = Tensioni;
    Fd360 = Fd;
    Fq360 = Fq;
    time = linspace(tempo(1),round(360/Cas.sim_angle)*tempo(end),length(Coppia360));
    theta360 = (Mac.th0 * pi/180 + time/1000 * Cas.n * pi/30 * Mac.p);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coppia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coppia = [Coppia sum(Ft.values(:,13:end),2)];
%Coppia  =   Ft.values(:,12:13);     %Nm

figure(1), hold on
% plot(tempo,abs(Coppia),'b-x')
% plot(tempo,abs(Coppia),'b-x')
% plot(tempo,Coppia - mean(Coppia),'b-x')
plot(theta360*180/pi,Coppia - mean(Coppia),'b-x')
hold off, grid on
legend('maxwell statore e rotore mediati')
xlabel('tempo [ms]'), ylabel('[Nm]')
title(['Coppia-' Descr ' - ' datestr(now)])
saveas(gcf,[FilePath 'fig1.fig'])

%%%%%%%%%%%%%%%%%%
%% Flussi
%%%%%%%%%%%%%%%%%%
figure(2), hold on
plot(tempo, Flussi(:,1),'b-x')
plot(tempo, Flussi(:,2),'r-x')
plot(tempo, Flussi(:,3),'g-x')
hold off, grid on
xlabel('tempo [ms]'), ylabel('[Vs]')
title(['Flusso-' Descr])
saveas(gcf,[FilePath 'fig2.fig'])



%%%%%%%%%%%%%%%%%%
%% EMF
%%%%%%%%%%%%%%%%%%
figure(3)
hold on
plot(tempo, Tensioni(:,1),'b-')
plot(tempo, Tensioni(:,2),'r-')
plot(tempo, Tensioni(:,3),'g-')
hold off
grid on
xlabel('tempo [ms]')
ylabel('[V]')
title(['EMF-' Descr])


%% %%%%%%%%%%%%%%%%%%%%%%%%%%
%     Stampa Elaborazione   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
hold on
plot(tempo,Fd,'b-',tempo,Fq,'r-x')
hold off
grid on
xlabel('tempo [ms]')
ylabel('[Vs]')
title(['Flussi D(blu) e Q(rosso)' Descr])
saveas(gcf,[FilePath 'fig4.fig'])



