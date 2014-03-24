
function geo=draw_rotor_3U(geo)
%%
%% INIZIALIZZAZIONE DATI DI INGRESSO:
xr = geo.xr;                    % Raggio del rotore al traferro
x0 = geo.x0;                    % Centro fittizio
rshaft = geo.Ar;                    % Raggio albero
Ar=geo.Ar;
l = geo.l;                      % Lunghezza pacco
g = geo.g;                      % Traferro
pont0 = geo.pont0;              % Ponticelli al traferro (i ponticelli al traferro hanno lo spessore di un arco lungo pont0)

p = geo.p;                      % Paia poli
nlay = geo.nlay;                % N° layers

alpha = geo.alpha;              % Angoli alpha (posizione barriere)
dalpha = geo.dalpha;            % Angoli dalpha
hc = geo.hc;                    % Altezze hc

racc_pont = geo.racc_pont;      % racc_pont=1*pont0 <- per i ponticelli radiali.
ang_pont0 = geo.ang_pont0;      % Ampiezza dell'angolo (in gradi) da spazzare con  raggio xr in modo da ottenre un arco lungo pont0

nmax = geo.nmax;                % Velocità max (rpm) per la valutazione della sollecitazione centrifuga più gravosa (-> ponticelli)

Delta_X=geo.Delta_X;            % Inizio parte rettilinea 1th layer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ampiezza di un arco lungo pont0
% ang_pont0 = pont0 / xr * 180/pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DISEGNO DELLE BARRIERE E DEI PONTICELLI
% INTRODUZIONE DI ALCUNE GRANDEZZE GEOMETRICHE UTILI E RIFERIMENTO AD UN CENTRO FITTIZIO DI COORDINATE (x0,0)

beta = 180/pi * calc_apertura_cerchio(pi/180*alpha,xr,x0);                  % La funzione calc_apertura_cerchio riceve in input le coordinate polari,
% (xr, alpha) (alpha in rad), di un punto generico. Queste sono calcolate
% rispetto al centro (0,0). In output la funzione restituisce l'apertura
% angolare (in rad) dello stesso punto rispetto al centro preso come
% riferimento (ha coordinate:(x0,0)).
% I punti di cui, in questo caso, si calcolano le aperture angolari rispetto
% al centro di riferimento sono i punti mediani delle barriere, presi in
% corrispondenza del traferro.

r = (x0 - xr * cos(alpha*pi/180))./(cos(beta*pi/180));                      % Di questi stessi punti, si calcolano anche le distanze dal centro (x0,0)
% e le si memorizzano nel vettore r.

r_all = [];                                                                 % Nel vettore r_all si memorizzano invece le distanze, sempre rispetto al
% centro di riferimento, dei punti, presi in corrispondenza del traferro,
for jj = 1:nlay                                                             % che individuano l'inizio e la fine delle nlay barriere.
    r_all = [r_all r(jj)-hc(jj)/2 r(jj)+hc(jj)/2];
end


%% Posizione banane su asse q (convenzione assi VAGATI)
% 2013/07/06 MG punto banana su asse q serve per i successivi export in DXF

% hcc=(r_all(2:2:end))-(r_all(1:2:end));
[xpont,ypont] = calc_intersezione_cerchi((xr-pont0), r, x0);
[xcbar,ycbar] = calc_intersezione_cerchi((xr-pont0-hc/2), r, x0);

hfe=[xr-xcbar(1)-hc(1)/2,r_all(3:2:end)-r_all(2:2:end-1),x0-r_all(end)]; %% 

rpont_x0=sqrt(ypont.^2+(x0-xpont).^2);
[alphapont,rpont] = cart2pol(xpont,ypont);
Bx0=x0-(rpont_x0);
%%
%% Calc di xpont, ypont:
for ii=1:nlay
[a,b,c]=retta_per_2pti(0,0,xcbar(ii),ycbar(ii));
A=1+b^2/a^2; B=2*b*c/a; C=(c^2/a^2-(xr-pont0)^2);
ytemp=roots([A,B,C]); ypont(ii)=ytemp(find(ytemp>=0));
xpont(ii)=-(b*ypont(ii)+c)/a;
end
% keyboard

for ii=1:nlay
    if ii==1
    B2k(ii)=xpont(ii)+Delta_X*((x0-r_all(1))-xpont(ii)); 
    B1k(ii)=B2k(ii)-hc(ii);
    
    else
    B1k(ii)=B1k(ii-1)-hfe(ii)-hc(ii); 
    B2k(ii)=B1k(ii-1)-hfe(ii);
    end

end
%% Controllo che si verifichi la configurazione I2U per Delta_X=1
if (Delta_X==0)
    I2Uconfig=1;
else
    I2Uconfig=0;
end

%% Intersezione circonferenze punti al traferro:
[xTraf2,yTraf2] = calc_intersezione_cerchi(xr-pont0, r_all(1:2:end), x0);
[xTraf1,yTraf1] = calc_intersezione_cerchi(xr-pont0, r_all(2:2:end), x0);

%% Intersezione tra rette che compongono i lati delle barriere di flux:
% retta verticale
for ii=1:nlay
    if (ii==1 && I2Uconfig==1)
    XpBar2(ii)=xpont(1);
    YpBar2(ii)=ypont(1);

    else
    a1=1;
    b1=0;
    c1=-B2k(ii);
    m2=tan(pi/2/p);
    a2=m2;
    b2=-1;
    c2=(yTraf2(ii)-m2*xTraf2(ii));
    [x y]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
    XpBar2(ii)=x;
    YpBar2(ii)=y;
    end
end

for ii=1:nlay
    if (ii==1 && I2Uconfig==1)
    XpBar1(ii)=B1k(1);
    YpBar1(ii)=ypont(1);

    else
    a1=1;
    b1=0;
    c1=-B1k(ii);
    m2=tan(pi/2/p);
    a2=m2;
    b2=-1;
    c2=(yTraf1(ii)-m2*xTraf1(ii));
    [x y]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
    XpBar1(ii)=x;
    YpBar1(ii)=y;
    end
end


for ii=1:nlay
    
    if (ii==1 && I2Uconfig==1)
        
    xD1k=B1k(1);
    yD1k=ycbar(1);
    xxD1k(ii)=B1k(1);
    yyD1k(ii)=ycbar(1);
    XpBar1(1)=xD1k;
    YpBar1(1)=yD1k;

    rcir=hc(1)/2;
    xcir=xcbar(1); ycir=ycbar(1);
    
    angT2=atan2(yD1k-ycir,xD1k-xcir);
    angT1=atan2(hc(1)/2,0);
    if (angT2<0)
        angT2=pi+angT2;
    end

    arcLayTraf1(ii)=(angT2-angT1)*180/pi;
    XcRibTraf1(ii)=xcbar(1);
    YcRibTraf1(ii)=ycbar(1);

    clear XCer YCer
%     XCer=xcir+rcir*cos(linspace(angT1,angT2,20));
%     YCer=ycir+rcir*sin(linspace(angT1,angT2,20));
%     figure(100);hold on;plot(XCer,YCer,'+m');hold off;
%     figure(100);hold on;plot(xD1k,yD1k,'rs');plot(xcbar(ii),ycbar(ii),'gs'); hold off;
%     figure(100);hold on;plot([0,xpont(ii)]',[0,ypont(ii)]','-b'); hold off;

    else
    [xD1k,yD1k,xc,yc,rc]=tg_cir(XpBar1(ii),YpBar1(ii),xTraf1(ii),yTraf1(ii),xpont(ii),ypont(ii));
    xxD1k(ii)=xD1k;
    yyD1k(ii)=yD1k;
    
    angT2=atan2(yD1k-yc,xD1k-xc);
    angT1=atan2(ypont(ii)-yc,xpont(ii)-xc);
    if (angT2<0)
        angT2=pi+angT2;
    end
    arcLayTraf1(ii)=(angT2-angT1)*180/pi;
    XcRibTraf1(ii)=xc;
    YcRibTraf1(ii)=yc;

    clear XCer YCer
%     XCer=xc+rc*cos(linspace(angT1,angT2,20));
%     YCer=yc+rc*sin(linspace(angT1,angT2,20));
%     figure(100);hold on;plot(XCer,YCer,'+m');hold off;
%     figure(100);hold on;plot(xD1k,yD1k,'rs');plot(xc,yc,'gs'); hold off;
%     figure(100);hold on;plot([0,xpont(ii)]',[0,ypont(ii)]','-b'); hold off;
    end
end

for ii=1:nlay
    if (ii==1&&I2Uconfig==1)
    xD2k=B2k(1);
    yD2k=ycbar(1);
    xxD2k(ii)=xD2k;
    yyD2k(ii)=yD2k;
    XpBar2(1)=xD2k;
    YpBar2(1)=yD2k;
    
    rcir=hc(1)/2;
    xcir=xcbar(1); ycir=ycbar(1);
    
    angT2=atan2(yD2k-ycir,xD2k-xcir);
    angT1=atan2(hc(1)/2,0);
    if (angT2<0)
        angT2=pi+angT2;
    end

    arcLayTraf2(ii)=(angT2-angT1)*180/pi;
    XcRibTraf2(ii)=xcbar(1);
    YcRibTraf2(ii)=ycbar(1);

    clear XCer YCer
%     XCer=xcir+rcir*cos(linspace(angT1,angT2,20));
%     YCer=ycir+rcir*sin(linspace(angT1,angT2,20));
%     figure(100);hold on;plot(XCer,YCer,'+m');hold off;
%     figure(100);hold on;plot(xD1k,yD1k,'rs');plot(xcbar(ii),ycbar(ii),'gs'); hold off;
%     figure(100);hold on;plot([0,xpont(ii)]',[0,ypont(ii)]','-b'); hold off;


    else
        
    [xD2k,yD2k,xc,yc,rc]=tg_cir(XpBar2(ii),YpBar2(ii),xTraf2(ii),yTraf2(ii),xpont(ii),ypont(ii));
    xxD2k(ii)=xD2k;
    yyD2k(ii)=yD2k;
    
    angT1=atan2(yD2k-yc,xD2k-xc);
    angT2=atan2(ypont(ii)-yc,xpont(ii)-xc);
    XcRibTraf2(ii)=xc;
    YcRibTraf2(ii)=yc;
    arcLayTraf2(ii)=(angT2-angT1)*180/pi;
    xxD2k(ii)=xD2k;
    yyD2k(ii)=yD2k;
        clear XCer YCer
%     XCer=xc+rc*cos(linspace(angT1,angT2,20));
%     YCer=yc+rc*sin(linspace(angT1,angT2,20));
%     figure(100);hold on;plot(XCer,YCer,'+m');hold off;
    end

end
%%
%% Controllo di sicurezza per evitare che la circonferenza del rib traf termini dopo il lato di barriera
%%
for ii=1:nlay
    if (xxD2k(ii)<=XpBar2(ii))
        if(XpBar2(ii)>xpont(ii))
            XpBar2(ii)=xpont(ii);
            B2k(ii)=xpont(ii);
        end
        rc=sqrt((xxD2k(ii)-XcRibTraf2(ii))^2+(yyD2k(ii)-YcRibTraf2(ii))^2);
        
        xxD2k(ii)=XpBar2(ii);
        yyD2k(ii)=-sqrt(rc^2-(XpBar2(ii)-XcRibTraf2(ii))^2)+YcRibTraf2(ii);
        angT1=atan2(yyD2k(ii)-YcRibTraf2(ii),xxD2k(ii)-XcRibTraf2(ii));
        angT2=atan2(ypont(ii)-YcRibTraf2(ii),xpont(ii)-XcRibTraf2(ii));
        arcLayTraf2(ii)=(angT2-angT1)*180/pi;
    end
end
% 

% %%
% xcir_plot=[0:0.5:xr];
% ycir_plot=sqrt(xr^2-xcir_plot.^2);
% figure(100);hold on;
% xo=x0-rpont_x0'*cos(0:0.1:pi/2);
% yo=rpont_x0'*sin(0:0.1:pi/2);
% plot(xo',yo','--r','LineWidth',2); axis([0 xr+0.5 0 xr+0.5]); axis square
% plot(B1k,0,'ob');plot(B2k,0,'ob');
% plot(xpont,ypont,'*c');
% plot(XpBar2,YpBar2,'bs');
% plot(XpBar1,YpBar1,'bs');
% plot(xTraf1,yTraf1,'*m');
% plot(xTraf2,yTraf2,'*m');
% plot(xD1k,yD1k,'r^');
% plot(xcir_plot,ycir_plot,'k');
% 
% for ii=1:nlay
%    plot([XpBar1(ii),xxD1k(ii)],[YpBar1(ii),yyD1k(ii)],'b','LineWidth',2); 
%    plot([XpBar2(ii),xxD2k(ii)],[YpBar2(ii),yyD2k(ii)],'b','LineWidth',2); 
%    plot([B1k(ii),XpBar1(ii)],[0,YpBar1(ii)],'b','LineWidth',2); 
%    plot([B2k(ii),XpBar2(ii)],[0,YpBar2(ii)],'b','LineWidth',2); 
% 
% end
% hold off;

%%
%%
%%   
% Determinazione dei punti caratteristici per il calcolo delle perdite nel
% ferro:
calcolo_posizioni_Pfe;
%% Valutazione ponticelli radiali:

calc_ribs_rad;
%% 
%% Salvataggio dei dati finali:

[xc,yc] = calc_intersezione_cerchi(xr-pont0-hc/2, r, x0);

geo.B1k=B1k;
geo.B2k=B2k;
geo.Bx0=Bx0;
geo.xpont=xpont;
geo.ypont=ypont;
geo.xxD1k=xxD1k;
geo.yyD1k=yyD1k;
geo.xxD2k=xxD2k;
geo.yyD2k=yyD2k;
geo.XpBar1=XpBar1;
geo.YpBar1=YpBar1;
geo.XpBar2=XpBar2;
geo.YpBar2=YpBar2;

geo.xTraf1=xTraf1;
geo.xTraf2=xTraf2;
geo.yTraf1=yTraf1;
geo.yTraf2=yTraf2;

geo.arcLayTraf1=arcLayTraf1;
geo.arcLayTraf2=arcLayTraf2;

geo.XcRibTraf1=XcRibTraf1;
geo.YcRibTraf1=YcRibTraf1;
geo.XcRibTraf2=XcRibTraf2;
geo.YcRibTraf2=YcRibTraf2;

geo.r_fe=rfe;
geo.x_fe=xFe;
geo.y_fe=yFe;
geo.xc=xc;
geo.yc=yc;

%% Punti per i ribs radiali
geo.XpontRadDx=XpontRadDx;
geo.YpontRadDx=YpontRadDx;
geo.XpontRadSx=XpontRadSx;
geo.YpontRadSx=YpontRadSx;
geo.XpontRadBarDx=XpontRadBarDx;
geo.XpontRadBarSx=XpontRadBarSx;
geo.YpontRadBarDx=YpontRadBarDx;
geo.YpontRadBarSx=YpontRadBarSx;

beta_f=0;
hf=[xr,B1k]-[B2k,Ar]; %calcolo dei Delta fi ferro di rotore
geo.hf = hf;
geo.beta_f = beta_f;
geo.pont=pont;
% keyboard
end