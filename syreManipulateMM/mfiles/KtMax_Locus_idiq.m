% 08 11 09

% I CONTOUR evaluation --> Kt max loci
% CurveIsoII is the I Levels Contour in the id-iq plane

cd, debug = 0;

Curve = CurveIsoII; Level = ILevel;
Tot = 0; Valori = 0; NuovoValore=[];

a = cell(length(Level),1);

id_KtMax = []; iq_KtMax = []; T_KtMax = [];

m = 1; finito = 0; count_finito = 0;

while ((Tot+Valori)<length(Curve) && finito == 0),
    
    Tot=Tot+Valori+1;
    Valori = Curve(2,Tot);    %n of points at this level
    
    if NuovoValore==Curve(1,Tot)
        Tot=Tot+Valori+1;
        Valori=Curve(2,Tot);
    end
    NuovoValore=Curve(1,Tot);
    
    idIso=Curve(1,Tot+[1:Valori]);
    iqIso=Curve(2,Tot+[1:Valori]);
    
    T_I = interp2(id,iq,TI,idIso,iqIso);    % T at given I - 1st round
    [Tmax_I, indiceM]=max(T_I);
    
    id_KtMax(m) = idIso(indiceM);           % luogo max T / I (Nm/A)
    iq_KtMax(m) = iqIso(indiceM);
    
    if ((indiceM == 1)||(strcmp(axes_type,'SR') && (indiceM == length(idIso))))
        count_finito = count_finito + 1;
        id_KtMax(m) = NaN;     % luogo max T / I (Nm/A)
        iq_KtMax(m) = NaN;
        Tmax_I = NaN;
    end
    
    if count_finito == 5
        finito = 1;
    end
    T_KtMax(m) = Tmax_I;
    
    %         figure(200)
    %         plot(idIso,iqIso,'-x'), axis equal, hold on
    %         plot(id_KtMax(m),iq_KtMax(m),'xr'), hold off
    %         title(['Iso corrente ' num2str(ILevel(m)) ' A - cerco la max coppia'])
    %         indiceM
    %         keyboard
    
    m = m + 1;
end

id_KtMax = id_KtMax(not(isnan(id_KtMax)));
iq_KtMax = iq_KtMax(not(isnan(iq_KtMax)));
T_KtMax = T_KtMax(not(isnan(T_KtMax)));

% poly fit KtMax
if strcmp(axes_type,'SR')
    %% IPM style axes
    [p_KtMax_i,s] = polyfit(id_KtMax,iq_KtMax,7);
    [p_KtMax_T,s] = polyfit(id_KtMax,T_KtMax,7);
    id_KtMax_p = linspace(0,max(id_KtMax)*1.05,length(id_KtMax));
    iq_KtMax_p = polyval(p_KtMax_i,id_KtMax_p);
    T_KtMax_p = polyval(p_KtMax_T,id_KtMax_p);
else
    %% SPM style axes
    [p_KtMax_i,s] = polyfit(iq_KtMax,id_KtMax,7);
    [p_KtMax_T,s] = polyfit(iq_KtMax,T_KtMax,7);
    iq_KtMax_p = linspace(0,max(iq_KtMax),length(iq_KtMax))
    id_KtMax_p = polyval(p_KtMax_i,iq_KtMax_p);
    T_KtMax_p = polyval(p_KtMax_T,iq_KtMax_p);
end

if not(exist('Imax'))
    Imax = max(max(II))/sqrt(2);
end
ind = find(abs(iq_KtMax_p + j* id_KtMax_p)>=0.94*Imax,1,'first');

if debug
    figure(3000)
    [c,h]=contour(id,iq,TI,linspace(0,T_KtMax_p(end),11)); clabel(c,h), hold on
    [c,h]=contour(id,iq,II,linspace(0,Imax,11)); clabel(c,h),
    % MTPA - punti
    plot(id_KtMax,iq_KtMax,'kx'),
    % MTPA - interpolata
    plot(id_KtMax_p,iq_KtMax_p,'b','LineWidth',[0.75])
    axis equal, hold off
    grid on, xlabel('Fd [Vs]'),ylabel('Fq [Vs]'), hold off
    title(['MTPA (KtMax) - Imax = ' num2str(Imax/sqrt(2),4) 'Arms'],'FontSize',[12],'FontWeight','bold');
    axis equal
end

%% uso le curve interpolate
if(0)
    id_KtMax = id_KtMax_p;
    iq_KtMax = iq_KtMax_p;
    T_KtMax = T_KtMax_p;
end
I_KtMax = abs(id_KtMax + j* iq_KtMax);
save([pathname 'ktMax_idiq'],'id_KtMax','iq_KtMax','T_KtMax')
if exist('dTpp')
    dTpp_KtMax = interp2(id,iq,dTpp,id_KtMax,iq_KtMax);    % dT pk pk on MTPA
    dTpp_KtMax = dTpp_KtMax(not(isnan(T_KtMax)));
    save ([pathname 'ktMax_idiq'],'dTpp_KtMax','-append');
end
% save('mat\ktMax_FdFq','fd_KtMax','fq_KtMax','T_KtMax')

% print(gcf,'-dpdf','-r300',['Mappe_IPMScooter\MotorMap_Flux\14_Ktmax_locus'])
