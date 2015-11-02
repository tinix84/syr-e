% 22 05 07

% I CONTOUR evaluation --> Kt max loci
% run by C_Map4Control_LdLq.m
% CurveIsoI is the I Levels Contour in the fd-fq plane

cd
Curve = CurveIsoI;
Level = ILevel;
Tot = 0;
Valori = 0;
NuovoValore=[];

a = cell(length(Level),1);

fd_KtMax = [];
fq_KtMax = [];
T_KtMax = [];

m = 1;
finito = 0;
count_finito = 0;

while (m<=length(ILevel) && finito == 0),

    Tot=Tot+Valori+1;
    Valori = Curve(2,Tot);    %n of points at this level

    if NuovoValore==Curve(1,Tot)
        Tot=Tot+Valori+1;
        Valori=Curve(2,Tot);
    end
    NuovoValore=Curve(1,Tot);

    fdIso=Curve(1,Tot+[1:Valori]);
    fqIso=Curve(2,Tot+[1:Valori]);

    T_I = interp2(FD,FQ,T,fdIso,fqIso);    % T at given I - 1st round
    [Tmax_I, indiceM]=max(T_I);

    %       fd1 = fdIso(indiceM-1);
    %       fd2 = fdIso(indiceM+1);
    %       d_fd = linspace(fd1,fd2,50);
    %       fq1 = fqIso(indiceM-1);
    %       fq2 = fqIso(indiceM+1);
    %       d_fq = linspace(fq1,fq2,50);
    %
    %       T_I = interp2(FD,FQ,T,d_fd,d_fq);    % T at given I - 2nd round
    %
    %       [Tmax_I, indiceM]=max(T_I);
    %       fd_KtMax(m) = d_fd(indiceM);     % luogo max T / I (Nm/A)
    %       fq_KtMax(m) = d_fq(indiceM);

    fd_KtMax(m) = fdIso(indiceM);     % luogo max T / I (Nm/A)
    fq_KtMax(m) = fqIso(indiceM);
    
    PF_KtMax(m) = interp2(fd,fq,PF,fd_KtMax(m),fq_KtMax(m));
    
    if (indiceM == 1 || indiceM == length(fdIso))
        count_finito = count_finito + 1;
    end
    if count_finito == 2
        finito = 1;
    end
    T_KtMax(m) = Tmax_I;

%     figure(200)
%     plot(fdIso,fqIso,'-x'), axis equal, hold on
%     plot(fd_KtMax(m),fq_KtMax(m),'xr'), hold off
%     title(['Iso corrente ' num2str(ILevel(m)) ' A - cerco la max coppia'])
%     indiceM
%     pause

    %deltaIso = atan(fqIso./fdIso)
    m = m + 1;
end

% poly fit KtMax
[p_KtMax_f,s] = polyfit(fd_KtMax,fq_KtMax,7);
[p_KtMax_T,s] = polyfit(fd_KtMax,T_KtMax,7);
[p_KtMax_PF,s] = polyfit(fd_KtMax,PF_KtMax,7);
fd_KtMax_p = linspace(0,1.05*max(fd_KtMax),length(fd_KtMax))
fq_KtMax_p = polyval(p_KtMax_f,fd_KtMax_p);
T_KtMax_p = polyval(p_KtMax_T,fd_KtMax_p);
PF_KtMax_p = polyval(p_KtMax_PF,fd_KtMax_p);

fig2 = figure;
temp = max(max(abs(F)));
[c,h]=contour(FD,FQ,F,temp * [0.1 0.5 0.9]);
clabel(c,h); hold on
temp = max(max(abs(T)));
[c,h]=contour(FD,FQ,T,temp * (0.1:0.1:0.9)); clabel(c,h);
[c,h]=contour(FD,FQ,I,Imax * [0.5 1]); clabel(c,h);
% MTPA - punti 
plot(fd_KtMax,fq_KtMax,'kx','LineWidth',[1.5])
% MTPA - interpolata
plot(fd_KtMax_p,fq_KtMax_p,'b','LineWidth',[0.75])
grid on, xlabel('Fd [Vs]'),ylabel('Fq [Vs]'), hold off
title(['MTPA (KtMax) - Imax = ' num2str(Imax/sqrt(2),4) 'Arms'],'FontSize',[12],'FontWeight','bold');
axis equal

% % debug
% figure(100)
% plot(fd_KtMax,T_KtMax,'kx'), hold on
% plot(fd_KtMax_p,T_KtMax_p), 
% plot(fd_KtMax,PF_KtMax,'kx'),
% plot(fd_KtMax_p,PF_KtMax_p), hold off, grid on


%% uso le curve interpolate
fd_KtMax = fd_KtMax_p;
fq_KtMax = fq_KtMax_p;
T_KtMax = T_KtMax_p;
PF_KtMax = PF_KtMax_p;

save([pathname 'ktMax_FdFq'],'fd_KtMax','fq_KtMax','T_KtMax','PF_KtMax')
% save('mat\ktMax_FdFq','fd_KtMax','fq_KtMax','T_KtMax')

% print(gcf,'-dpdf','-r300',['Mappe_IPMScooter\MotorMap_Flux\14_Ktmax_locus'])
