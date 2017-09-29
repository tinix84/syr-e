% 01 01 2012

% F CONTOUR evaluation --> Kv max


% cd
Curve = CurveIsoTI;
Level = TLevel;
Tot = 0;
Valori = 0;
NuovoValore=[];

a = cell(length(Level),1);

id_KvMax = [];
iq_KvMax = [];
F_KvMax = [];

m = 1;
finito = 0;
count_finito = 0;

while (m<=length(Level) && finito == 0),
    
    Tot=Tot+Valori+1;
    Valori = Curve(2,Tot);    %n of points at this level
    
    if NuovoValore==Curve(1,Tot)
        Tot=Tot+Valori+1;
        Valori=Curve(2,Tot);
    end
    NuovoValore=Curve(1,Tot);
    
    idIso=Curve(1,Tot+[1:Valori]);
    iqIso=Curve(2,Tot+[1:Valori]);
    
    fdIso = interp2(Id,Iq,Fd,idIso,iqIso);
    fqIso = interp2(Id,Iq,Fq,idIso,iqIso);
    
    [FtangT, indiceM]=min(fdIso.^2 + fqIso.^2);
    id_KvMax(m) = idIso(indiceM);          % luogo max T / Flux (Nm/Vs)
    iq_KvMax(m) = iqIso(indiceM);
    fd_KvMax(m) = fdIso(indiceM);          % luogo max T / Flux (Nm/Vs)
    fq_KvMax(m) = fqIso(indiceM);
    
    F_KvMax(m) = sqrt(fdIso(indiceM)^2 + fqIso(indiceM)^2);
    
    T_KvMax(m) = Level(m);
    
    %         figure(200)
    %         plot(idIso,iqIso,'-x'), axis equal, hold on
    %         plot(id_KvMax(m),iq_KvMax(m),'xr'), hold off
    %         title(['Iso corrente ' num2str(ILevel(m)) ' A - cerco la max coppia'])
    %         indiceM
    %         pause
    %
    m = m + 1;
end

%% elimina i bordi della MTPV (fuori dal dominio di flusso)
test = id_KvMax == max(id);
id_KvMax = id_KvMax(test == 0);
iq_KvMax = iq_KvMax(test == 0);
T_KvMax = T_KvMax(test == 0);
fd_KvMax = fd_KvMax(test == 0);
fq_KvMax = fq_KvMax(test == 0);
F_KvMax = F_KvMax(test == 0);

test = id_KvMax == min(id);
id_KvMax = id_KvMax(test == 0);
iq_KvMax = iq_KvMax(test == 0);
T_KvMax = T_KvMax(test == 0);
fd_KvMax = fd_KvMax(test == 0);
fq_KvMax = fq_KvMax(test == 0);
F_KvMax = F_KvMax(test == 0);

test = iq_KvMax == max(iq);
id_KvMax = id_KvMax(test == 0);
iq_KvMax = iq_KvMax(test == 0);
T_KvMax = T_KvMax(test == 0);
fd_KvMax = fd_KvMax(test == 0);
fq_KvMax = fq_KvMax(test == 0);
F_KvMax = F_KvMax(test == 0);

test = iq_KvMax == min(iq);
id_KvMax = id_KvMax(test == 0);
iq_KvMax = iq_KvMax(test == 0);
T_KvMax = T_KvMax(test == 0);
fd_KvMax = fd_KvMax(test == 0);
fq_KvMax = fq_KvMax(test == 0);
F_KvMax = F_KvMax(test == 0);

% no PM machines
if fm == 0
    id_KvMax = [0 id_KvMax];
    iq_KvMax = [0 iq_KvMax];
    T_KvMax = [0 T_KvMax];
    fd_KvMax = fd_KvMax(test == 0);
    fq_KvMax = fq_KvMax(test == 0);
    F_KvMax = F_KvMax(test == 0);
end

if not(isempty(id_KvMax))
    % poly fit KvMax
    if strcmp(axes_type,'SR')
        %% SR style axes
        [p_KvMax_i,s] = polyfit(id_KvMax,iq_KvMax,7);
        [p_KvMax_T,s] = polyfit(id_KvMax,T_KvMax,7);
        [p_KvMax_F,s] = polyfit(id_KvMax,F_KvMax,7);
        [p_KvMax_fd,s] = polyfit(id_KvMax,fd_KvMax,7);
        [p_KvMax_fq,s] = polyfit(id_KvMax,fq_KvMax,7);
        id_KvMax_p = linspace(0,max(id_KvMax)*1.05,length(id_KvMax));
        iq_KvMax_p = polyval(p_KvMax_i,id_KvMax_p);
        T_KvMax_p = polyval(p_KvMax_T,id_KvMax_p);
        F_KvMax_p = polyval(p_KvMax_F,id_KvMax_p);
        fd_KvMax_p = polyval(p_KvMax_fd,id_KvMax_p);
        fq_KvMax_p = polyval(p_KvMax_fq,id_KvMax_p);
    else
        %% PM style axes
        [p_KvMax_i,s] = polyfit(iq_KvMax,id_KvMax,7);
        [p_KvMax_T,s] = polyfit(iq_KvMax,T_KvMax,7);
        [p_KvMax_F,s] = polyfit(id_KvMax,F_KvMax,7);
        [p_KvMax_fd,s] = polyfit(id_KvMax,fd_KvMax,7);
        [p_KvMax_fq,s] = polyfit(id_KvMax,fq_KvMax,7)
        iq_KvMax_p = linspace(0,max(iq_KvMax),length(iq_KvMax));
        id_KvMax_p = polyval(p_KvMax_i,iq_KvMax_p);
        T_KvMax_p = polyval(p_KvMax_T,iq_KvMax_p);
        F_KvMax_p = polyval(p_KvMax_F,id_KvMax_p);
        fd_KvMax_p = polyval(p_KvMax_fd,id_KvMax_p);
        fq_KvMax_p = polyval(p_KvMax_fq,id_KvMax_p);
    end
    
    if debug
        figure(4000)
        [c,h]=contour(id,iq,TI,max(max(TI))*[-0.5:0.1:0.9]); clabel(c,h), hold on
        [c,h]=contour(id,iq,FI,max(max(FI))*(0.1:0.1:0.9)); clabel(c,h),
        % MTPV - punti + interpolata
        plot(id_KvMax,iq_KvMax,'kx'), plot(id_KvMax_p,iq_KvMax_p,'b','LineWidth',[0.75])
        % MTPA - punti + interpolata
        plot(id_KtMax,iq_KtMax,'kx'), plot(id_KtMax_p,iq_KtMax_p,'b','LineWidth',[0.75])
        
        axis equal, hold off
        grid on, xlabel('id [A]'),ylabel('iq [A]'), hold off
        title(['MTPV (KvMax) - Imax = ' num2str(Imax/sqrt(2),4) 'Arms'],'FontSize',[12],'FontWeight','bold');
        axis equal
    end
    
    %% uso le curve interpolate (solo per IPM)
    if(0)
        id_KvMax = id_KvMax_p;
        iq_KvMax = iq_KvMax_p;
        T_KvMax = T_KvMax_p;
        I_KvMax = abs(id_KvMax + j* iq_KvMax);
        fd_KvMax = fd_KvMax_p;
        fq_KvMax = fq_KvMax_p;
        F_KvMax = F_KvMax_p;
    end
    
end

% PF_KvMax = PF_KvMax_p;
save([pathname 'kvMax_idiq'],'id_KvMax','iq_KvMax','F_KvMax','T_KvMax')
save([pathname 'kvMax_FdFq'],'fd_KvMax','fq_KvMax','F_KvMax','T_KvMax')

% save('mat\ktMax_FdFq','fd_KvMax','fq_KvMax','T_KvMax')

% print(gcf,'-dpdf','-r300',['Mappe_IPMScooter\MotorMap_Flux\14_Ktmax_locus'])
