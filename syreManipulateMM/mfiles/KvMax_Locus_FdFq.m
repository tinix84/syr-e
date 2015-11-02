% 25 05 07

% Kv max locus
% run by C_Map4Control_LdLq.m
% luogo b = 0 corrisponde a tangenza tra isoFlux e isoTorque

Curve = CurveIsoT;
Level = TLevel;
Tot = 0;
Valori = 0;
NuovoValore=[];

a = cell(length(Level),1);

fd_KvMax = [];
fq_KvMax = [];
F_KvMax = [];

fd_PFMax = [];
fq_PFMax = [];
PF_Max = [];

m = 1;
finito = 0;
count_finito = 0;

% while (m<=length(TLevel) && finito == 0),
while ((Tot+Valori)<length(Curve) && finito == 0),
    
    Tot=Tot+Valori+1;
    Valori = Curve(2,Tot);    %n of points at this level
    
    if NuovoValore==Curve(1,Tot)
        Tot=Tot+Valori+1;
        Valori=Curve(2,Tot);
    end
    NuovoValore=Curve(1,Tot);
    
    fdIso=Curve(1,Tot+[1:Valori]);
    fqIso=Curve(2,Tot+[1:Valori]);
    
    [FtangT, indiceM]=min(fdIso.^2 + fqIso.^2);
    fd_KvMax(m) = fdIso(indiceM);          % luogo max T / Flux (Nm/Vs)
    fq_KvMax(m) = fqIso(indiceM);
    
    PF_ = interp2(fd,fq,PF,fdIso,fqIso);
    [PF_max, indicePF_max]=max(PF_);
    
    fd_PFMax(m) = fdIso(indicePF_max);     % luogo max PF / T (PF/Nm)
    fq_PFMax(m) = fqIso(indicePF_max);
    
    PF_KvMax(m) = interp2(fd,fq,PF,fd_KvMax(m),fq_KvMax(m));
    
    F_KvMax(m) = sqrt(fdIso(indiceM)^2 + fqIso(indiceM)^2);
    PF_Max(m) = PF_max;
    
    %     deltaIso = atan(fqIso./fdIso)
    %     figure(200)
    %     plot(fdIso,fqIso), grid on, axis equal, hold on
    %     plot(fd_KvMax(m),fq_KvMax(m),'xr'),
    
    %     if m > 0
    %     indiceM
    %     length(fdIso)
    %     keyboard
    %     end
    
    %     if ((indiceM == 1))
    %         count_finito = count_finito + 1;
    %     end
    %
    %     if count_finito == 2
    %         finito = 1;
    %     end
    m = m +1;
    
    %     if (indiceM == 1 || indiceM == length(fdIso))
    %         finito = 1;
    %         %% DELETE LAST POINT
    %         fd_KvMax = fd_KvMax(1:end-1);
    %         fq_KvMax = fq_KvMax(1:end-1);
    %         F_KvMax = F_KvMax(1:end-1);
    %         PF_KvMax = PF_KvMax(1:end-1);
    %         PF_Max = PF_Max(1:end-1);
    %     end
    
end

T_KvMax = TLevel(1:m-1);
T_PFMax = TLevel(1:m-1);

if finito
    T_KvMax = TLevel(1:m-2);
    T_PFMax = TLevel(1:m-2);
end

%% elimina i bordi della MTPV (fuori dal dominio di flusso)
test = fd_KvMax == max(fd);
fd_KvMax = fd_KvMax(test == 0);
fq_KvMax = fq_KvMax(test == 0);
T_KvMax = T_KvMax(test == 0);

test = fd_KvMax == min(fd);
fd_KvMax = fd_KvMax(test == 0);
fq_KvMax = fq_KvMax(test == 0);
T_KvMax = T_KvMax(test == 0);

test = fq_KvMax == max(fq);
fd_KvMax = fd_KvMax(test == 0);
fq_KvMax = fq_KvMax(test == 0);
T_KvMax = T_KvMax(test == 0);

test = fq_KvMax == min(fq);
fd_KvMax = fd_KvMax(test == 0);
fq_KvMax = fq_KvMax(test == 0);
T_KvMax = T_KvMax(test == 0);

if debug
    fig2 = figure;
    [c,h]=contour(FD,FQ,F,max(max(F))*(0.1:0.1:0.9)); clabel(c,h); hold on
    [c,h]=contour(FD,FQ,T,max(max(TI))*(0.1:0.1:0.9));
    clabel(c,h);
    [c,h]=contour(FD,FQ,I,Imax * [0.5 1]);
    plot(fd_KvMax,fq_KvMax,'kx','LineWidth',[1.5]);
    grid on, xlabel('Fd [Vs]'),ylabel('Fq [Vs]')
    title(['MTPV (KvMax) - Imax = ' num2str(Imax/sqrt(2),4) 'Arms'],'FontSize',[12],'FontWeight','bold');
    % xlim([0 0.5]);ylim([-0.25 0.25])
    axis equal
end

% save('Mat\kvMax_FdFq','fd_KvMax','fq_KvMax','F_KvMax','T_KvMax')
save([pathname 'kvMax_FdFq'],'fd_KvMax','fq_KvMax','F_KvMax','T_KvMax','PF_KvMax')
save([pathname 'kPFMax_FdFq'],'fd_PFMax','fq_PFMax','PF_Max','T_PFMax')
%
% print(gcf,'-dpdf','-r300',['Mappe_IPMScooter\15_Kvmax_locus'])
