function [geo,temp]=nodes_rotor_Seg(geo)
%%
%% INIZIALIZZAZIONE DATI DI INGRESSO:
r = geo.r;                    % Raggio del rotore al traferro
x0 = geo.x0;                    % Centro fittizio
rshaft = geo.Ar;                    % Raggio albero
Ar=geo.Ar;
l = geo.l;                      % Lunghezza pacco
g = geo.g;                      % Traferro
pont0 = geo.pont0;              % Ponticelli al traferro (i ponticelli al traferro hanno lo spessore di un arco lungo pont0)

p = geo.p;                      % Paia poli
nlay = geo.nlay;                % N° layers

dalpha = geo.dalpha;            % Angoli dalpha
% Eval alpha
alpha = integral_fx(1:length(dalpha),dalpha);
geo.alpha=alpha;
dx=geo.dx;

racc_pont = geo.racc_pont;      % racc_pont=1*pont0 <- per i ponticelli radiali.
ang_pont0 = geo.ang_pont0;      % Ampiezza dell'angolo (in gradi) da spazzare con  raggio r in modo da ottenre un arco lungo pont0

nmax = geo.nmax;                % Velocità max (rpm) per la valutazione della sollecitazione centrifuga più gravosa (-> ponticelli)
sigma_max=geo.sigma_max;
hfe_min=geo.hfe_min;

alphaRot=alpha*pi/180+pi/2/p;   % angoli barriere di flux al traferro meccanici

rCirRib=1*pont0;
rlim=r;
DTrasl=geo.DTrasl;
% DTrasl=0;
%% Possible inital definition of
arcLayTraf1=zeros(1,nlay);
arcLayTraf2=zeros(1,nlay);
XcRibTraf1=zeros(1,nlay);
XcRibTraf2=zeros(1,nlay);
YcRibTraf1=zeros(1,nlay);
YcRibTraf2=zeros(1,nlay);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ampiezza di un arco lungo pont0
% ang_pont0 = pont0 / r * 180/pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DISEGNO DELLE BARRIERE E DEI PONTICELLI
% INTRODUZIONE DI ALCUNE GRANDEZZE GEOMETRICHE UTILI E RIFERIMENTO AD UN CENTRO FITTIZIO DI COORDINATE (x0,0)

beta = 180/pi * calc_apertura_cerchio(pi/180*alpha,r,x0);                  % La funzione calc_apertura_cerchio riceve in input le coordinate polari,
% (r, alpha) (alpha in rad), di un punto generico. Queste sono calcolate
% rispetto al centro (0,0). In output la funzione restituisce l'apertura
% angolare (in rad) dello stesso punto rispetto al centro preso come
% riferimento (ha coordinate:(x0,0)).
% I punti di cui, in questo caso, si calcolano le aperture angolari rispetto
% al centro di riferimento sono i punti mediani delle barriere, presi in
% corrispondenza del traferro.

rbeta = (x0 - r * cos(alpha*pi/180))./(cos(beta*pi/180));                      % Di questi stessi punti, si calcolano anche le distanze dal centro (x0,0)
% e le si memorizzano nel vettore rbeta.

[xpont,ypont] = calc_intersezione_cerchi(r-pont0, rbeta, x0);

LowDimBarrier=zeros(1,nlay);
for ii=1:nlay
    if (not(isreal(xpont(ii)))||not(isreal(ypont(ii))))
        %         betaTmp=atan(-((r-pont0)*sin(alpha(ii)*pi/180)/((r-pont0)*cos(alpha(ii)*pi/180)-x0)));
        xpont(ii)=(r-2*pont0)*cos(alpha(ii)*pi/180);
        ypont(ii)=(r-2*pont0)*sin(alpha(ii)*pi/180);
        LowDimBarrier(ii)=1;
    end
end

rpont_x0=sqrt(ypont.^2+(x0-xpont).^2);
[alphapont,rpont] = cart2pol(xpont,ypont);
Bx0=x0-(rpont_x0); geo.Bx0=Bx0;
%% Determination of air thickness and check the feasibility of the geometry
geo.Bx0=Bx0; % Initialization of central non-moved line of the flux barrier
geo = calcHcCheckGeoControlwDx(geo);
B1k=geo.B1k; B2k=geo.B2k;

[xc_temp,yc_temp]=calc_intersezione_cerchi(r,rbeta(nlay),x0);
dPointEndBar=calc_distanza_punti([xc_temp,yc_temp],[r*cos(pi/2/p), r*sin(pi/2/p)]);
if (dPointEndBar<(Bx0(nlay)-B1k(nlay)))
    B1k(nlay)=Bx0(nlay)-dPointEndBar+hfe_min/2;
end

hc=B2k-B1k;

ptmp=find(abs(hc)<pont0/2);
B2k(ptmp)=min([B1k(ptmp),B2k(ptmp)])+pont0/2;
clear ptmp;
%% Intersezione circonferenze punti al traferro:
[xTraf2,yTraf2] = calc_intersezione_cerchi(r-pont0, x0-B2k, x0);
[xTraf1,yTraf1] = calc_intersezione_cerchi(r-pont0, x0-B1k, x0);

for ii=1:nlay
    if (not(isreal(xTraf1(ii)))||not(isreal(yTraf1(ii)))||not(isreal(xTraf2(ii)))||not(isreal(yTraf2(ii)))||hc(ii)<0)
        %         betaTmp=atan(-((r-pont0)*sin(alpha(ii)*pi/180)/((r-pont0)*cos(alpha(ii)*pi/180)-x0)));
        xpont(ii)=(r-2*pont0)*cos(alpha(ii)*pi/180);
        ypont(ii)=(r-2*pont0)*sin(alpha(ii)*pi/180);
        LowDimBarrier(ii)=1;
    end
end
% keyboard
% Ktrasl=0.5;
% DXtrasl=Ktrasl*(xpont(1)-geo.B2k(1));
% B1k=B1k+DXtrasl;
% B2k=B2k+DXtrasl;

%% 1° barriera di flux: la prima barrira di flusso presenta raccordo
%% tondeggiante al traferro diviso in parti uguali pari alla più piccola
%% delle distanze da pont0.
dhoriz=calc_distanza_punti([B2k(1),0],[xTraf2(1),0]);
dDiag=calc_distanza_punti([B2k(1),0],[xTraf2(1),yTraf2(1)]);
if (dhoriz<=dDiag/4)
    xtraf1=xTraf1(1);
    ytraf1=yTraf1(1);
    xtraf2=xTraf2(1);
    ytraf2=yTraf2(1);
    
    dt11=abs((xpont(1)-xtraf1)+1j*(ytraf1-ypont(1)));
    dt12=abs((xpont(1)-xtraf2)+1j*(ytraf2-ypont(1)));
    [dt1,posMin]=min([dt11,dt12]);
    m=tan(pi/2/p);
    B=2*(ypont(1)-m*xpont(1));
    C=ypont(1)^2+m^2*xpont(1)^2-2*m*xpont(1)*ypont(1)-dt1^2*(1+m^2);
    q=B/2-sqrt(B^2/4-C);
    mprimo=-1/m;
    qprimo=ypont(1)-mprimo*xpont(1);
    if (posMin==1)
        xTraf2(1)=(qprimo-q)/(m-mprimo);
        yTraf2(1)=m*xTraf2(1)+q;
        xTraf1(1)=xtraf1(1);
        yTraf1(1)=ytraf1(1);
    else
        xTraf2(1)=(qprimo-q)/(m-mprimo);
        yTraf2(1)=m*xTraf2(1)+q;
        xTraf1(1)=xtraf1(1);
        yTraf1(1)=ytraf1(1);
        
    end
end
%% Control of the distance between (xTraf1,2;yTraf1,2) to mantein a minimum iron acceptable...
for ii=1:nlay-1
    dPto=calc_distanza_punti([xTraf2(ii+1),yTraf2(ii+1)],[xTraf1(ii),yTraf1(ii)]);
    dP1=calc_distanza_punti([xTraf1(ii),yTraf1(ii)],[xpont(ii),ypont(ii)]);
    dP2=calc_distanza_punti([xTraf2(ii+1),yTraf2(ii+1)],[xpont(ii+1),ypont(ii+1)]);
    if(dPto<0.0*hfe_min)
        if (dP1>dP2)
            [a,b,c]=retta_per_2pti(xTraf2(ii+1),yTraf2(ii+1),xTraf1(ii),yTraf1(ii));
            m=-a/b;
            A=1;
            B=-2*xTraf2(ii+1);
            C=xTraf2(ii+1)^2-(1.5*hfe_min)^2/(m^2+1);
            x_temp=roots([A,B,C]);
            xTraf1(ii)=x_temp(1);
            yTraf1(ii)=m*xTraf1(ii)-c/b;
        else
            [a,b,c]=retta_per_2pti(xTraf2(ii+1),yTraf2(ii+1),xTraf1(ii),yTraf1(ii));
            m=-a/b;
            A=1;
            B=-2*xTraf1(ii);
            C=xTraf1(ii)^2-(1.5*hfe_min)^2/(m^2+1);
            x_temp=roots([A,B,C]);
            xTraf2(ii+1)=x_temp(2);
            yTraf2(ii+1)=m*xTraf2(ii+1)-c/b;
        end
    end
end
%%
%% Intersezione tra rette che compongono i lati delle barriere di flux:
%%
% retta verticale
for ii=1:nlay
    if LowDimBarrier(ii)==1
        XpBar2(ii)=B2k(ii);
        YpBar2(ii)=0;
    else
        a1=1;
        b1=0;
        c1=-B2k(ii);
        m2=tan(pi/2/p);
        a2=m2;
        b2=-1;
        c2=(yTraf2(ii)-m2*xTraf2(ii));
        [x,y]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
        XpBar2(ii)=x;
        YpBar2(ii)=y;
    end
end

for ii=1:nlay
    if LowDimBarrier(ii)==1
        XpBar1(ii)=B1k(ii);
        YpBar1(ii)=0;
        
    else
        a1=1;
        b1=0;
        c1=-B1k(ii);
        m2=tan(pi/2/p);
        a2=m2;
        b2=-1;
        c2=(yTraf1(ii)-m2*xTraf1(ii));
        [x,y]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
        XpBar1(ii)=x;
        YpBar1(ii)=y;
    end
end

for ii=1:nlay
    if LowDimBarrier(ii)==1
        xxD1k(ii)=B1k(ii);
        yyD1k(ii)=0;
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
        XCer=xc+rc*cos(linspace(angT1,angT2,20));
        YCer=yc+rc*sin(linspace(angT1,angT2,20));
        %         figure(100);hold on;plot(XCer,YCer,'+m');hold off;
        %         figure(100);hold on;plot(xD1k,yD1k,'rs');plot(xc,yc,'gs'); hold off;
        %         figure(100);hold on;plot([0,xpont(ii)]',[0,ypont(ii)]','-b'); hold off;
    end
end

for ii=1:nlay
    if LowDimBarrier(ii)==1
        xxD2k(ii)=B2k(ii);
        yyD2k(ii)=0;
        
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
        XCer=xc+rc*cos(linspace(angT1,angT2,20));
        YCer=yc+rc*sin(linspace(angT1,angT2,20));
        %         figure(100);hold on;plot(XCer,YCer,'+m');hold off;
    end
end
%% Controllo di sicurezza per evitare che la circonferenza del rib traf termini dopo il lato di barriera

for ii=1:nlay
    if LowDimBarrier(ii)==0
        if (xxD2k(ii)<=XpBar2(ii))
            if(XpBar2(ii)>xpont(ii))
                XpBar2(ii)=xpont(ii);
                B2k(ii)=xpont(ii);
            end
            rc=sqrt((xxD2k(ii)-XcRibTraf2(ii))^2+(yyD2k(ii)-YcRibTraf2(ii))^2);
            
            xxD2k(ii)=XpBar2(ii);
            yyD2k(ii)=-sqrt(rc^2-(XpBar2(ii)-XcRibTraf2(ii))^2)+YcRibTraf2(ii);
            YpBar2(ii)=yyD2k(ii);
            angT1=atan2(yyD2k(ii)-YcRibTraf2(ii),xxD2k(ii)-XcRibTraf2(ii));
            angT2=atan2(ypont(ii)-YcRibTraf2(ii),xpont(ii)-XcRibTraf2(ii));
            arcLayTraf2(ii)=(angT2-angT1)*180/pi;
        end
    end
end
for ii=1:nlay
    if LowDimBarrier(ii)==0
        if (xxD1k(ii)<=XpBar1(ii))
            
            rc=sqrt((xxD1k(ii)-XcRibTraf1(ii))^2+(yyD1k(ii)-YcRibTraf1(ii))^2);
            
            xxD1k(ii)=XpBar1(ii);
            yyD1k(ii)=sqrt(rc^2-(XpBar1(ii)-XcRibTraf1(ii))^2)+YcRibTraf1(ii);
            YpBar1(ii)=yyD1k(ii);
            angT1=atan2(YcRibTraf1(ii)-yyD1k(ii),xxD1k(ii)-XcRibTraf1(ii));
            angT2=atan2(ypont(ii)-YcRibTraf1(ii),xpont(ii)-XcRibTraf1(ii));
            arcLayTraf1(ii)=(angT2-angT1)*180/pi;
        end
    end
end

if(xxD2k(1)<xpont(1) && yyD2k(1)>=ypont(1))
    YpBar2(1)=YpBar2(1)/2;
    yyD2k(1)=YpBar2(1);
end
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Raccordo barriere di flux:
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Si esegue un raccordo circolare al posto dello spigolo
XcRacc_B1=NaN*ones(1,nlay);
YcRacc_B1=NaN*ones(1,nlay);
xRaccR1_B1=NaN*ones(1,nlay);
yRaccR1_B1=NaN*ones(1,nlay);
xRaccR2_B1=NaN*ones(1,nlay);
yRaccR2_B1=NaN*ones(1,nlay);

for ii=1:length(B1k)
    
    if(XpBar1(ii)<xxD1k(ii) && strcmp(geo.RaccBarrier,'ON'))
        % Calcolo approssimativo della lunghezza media di barriera:
        d1B1 = calc_distanza_punti([B1k(ii),zeros(1,nlay)],[XpBar1(ii),YpBar1(ii)]);
        d2B1 = calc_distanza_punti([xxD1k(ii),yyD1k(ii)],[XpBar1(ii),YpBar1(ii)]);
        d1B2 = calc_distanza_punti([B2k(ii),zeros(1,nlay)],[XpBar2(ii),YpBar2(ii)]);
        d2B2 = calc_distanza_punti([xxD2k(ii),yyD2k(ii)],[XpBar2(ii),YpBar2(ii)]);
        dB1=sum([d1B1,d2B1]); dB2=sum([d1B2,d2B2]); d=mean([dB1,dB2]);
        % Il raggio è variabile per le varie barriere e corrisponde ad 1/4
        % della lunghezza media di barriera...
        R_RaccB1(ii)=d/4;
        [x,y,xtan1,ytan1,xtan2,ytan2]=cir_tg_2rette(B1k(ii),0,XpBar1(ii),YpBar1(ii),XpBar1(ii),YpBar1(ii),xxD1k(ii),yyD1k(ii),R_RaccB1(ii));
        XcRacc_B1(ii)=x(1); YcRacc_B1(ii)=y(1);
        xRaccR1_B1(ii)=xtan1(1);
        yRaccR1_B1(ii)=ytan1(1);
        xRaccR2_B1(ii)=xtan2(1);
        yRaccR2_B1(ii)=ytan2(1);
        
        %% Correct drawing control
        %% the following script control the feasibility of the arc, if not the radius is reduced to a feasible value
        if (xRaccR2_B1(ii)>=xxD1k(ii))
            d = calc_distanza_punti([xxD1k(ii),yyD1k(ii)],[XpBar1(ii),YpBar1(ii)]);
            R_RaccB1(ii)=d/2;
            [x,y,xtan1,ytan1,xtan2,ytan2]=cir_tg_2rette(B1k(ii),0,XpBar1(ii),YpBar1(ii),XpBar1(ii),YpBar1(ii),xxD1k(ii),yyD1k(ii),R_RaccB1(ii));
            XcRacc_B1(ii)=x(1); YcRacc_B1(ii)=y(1);
            xRaccR1_B1(ii)=xtan1(1);
            yRaccR1_B1(ii)=ytan1(1);
            xRaccR2_B1(ii)=xtan2(1);
            yRaccR2_B1(ii)=ytan2(1);
            
        end
        %%
        angT2=atan2(yRaccR2_B1(ii)-YcRacc_B1(ii),xRaccR2_B1(ii)-XcRacc_B1(ii));
        angT1=atan2(yRaccR1_B1(ii)-YcRacc_B1(ii),XcRacc_B1(ii)-xRaccR1_B1(ii))+3.14;
        if (angT2<0)
            angT2=pi+angT2;
        end
        
        clear XCer YCer
        XCer=XcRacc_B1(ii)+R_RaccB1(ii)*cos(linspace(angT1,angT2,20));
        YCer=YcRacc_B1(ii)+R_RaccB1(ii)*sin(linspace(angT1,angT2,20));
        %     figure(100);hold on;plot(XCer,YCer,'+m');hold off;
        %     figure(100);hold on;plot(xRaccR1_B1(ii),yRaccR1_B1(ii),'ms');hold off;
        %     figure(100);hold on;plot(xRaccR2_B1(ii),yRaccR2_B1(ii),'ms');hold off;
    else
        R_RaccB1(ii)=0;
    end
end
%%
XcRacc_B2=NaN*ones(1,nlay);
YcRacc_B2=NaN*ones(1,nlay);
xRaccR1_B2=NaN*ones(1,nlay);
yRaccR1_B2=NaN*ones(1,nlay);
xRaccR2_B2=NaN*ones(1,nlay);
yRaccR2_B2=NaN*ones(1,nlay);

for ii=1:length(B2k)
    
    if(XpBar2(ii)< xxD2k(ii) && strcmp(geo.RaccBarrier,'ON'))
        R_RaccB2(ii)=R_RaccB1(ii);
        [x,y,xtan1,ytan1,xtan2,ytan2]=cir_tg_2rette(B2k(ii),0,XpBar2(ii),YpBar2(ii),XpBar2(ii),YpBar2(ii),xxD2k(ii),yyD2k(ii),R_RaccB2(ii));
        % This control is necessary becouse in some cases there is an
        % intersection of the boundary of flux barrier due to the junctions,
        % statisitcally this seems to happen for the first flux barrier...
        if (ii>1 && abs((x(1)-XpBar1(ii-1))+1j*(y(1)-YpBar1(ii-1)))>R_RaccB2(ii))
            R_RaccB2(ii)=abs((x(1)-XpBar1(ii-1))+1j*(y(1)-YpBar1(ii-1)))+hfe_min;
            [x,y,xtan1,ytan1,xtan2,ytan2]=cir_tg_2rette(B2k(ii),0,XpBar2(ii),YpBar2(ii),XpBar2(ii),YpBar2(ii),xxD2k(ii),yyD2k(ii),R_RaccB2(ii));
        end
        
        XcRacc_B2(ii)=x(1); YcRacc_B2(ii)=y(1);
        xRaccR1_B2(ii)=xtan1(1);
        yRaccR1_B2(ii)=ytan1(1);
        xRaccR2_B2(ii)=xtan2(1);
        yRaccR2_B2(ii)=ytan2(1);
        
        %% Correct drawing control
        %% the following script control the feasibility of the arc, if not the radius is reduced to a feasible value
        if (xRaccR2_B2(ii)>=xxD2k(ii))
            d = calc_distanza_punti([xxD2k(ii),yyD2k(ii)],[XpBar2(ii),YpBar2(ii)]);
            R_RaccB2(ii)=d/4;
            [x,y,xtan1,ytan1,xtan2,ytan2]=cir_tg_2rette(B2k(ii),0,XpBar2(ii),YpBar2(ii),XpBar2(ii),YpBar2(ii),xxD2k(ii),yyD2k(ii),R_RaccB2(ii));
            XcRacc_B2(ii)=x(1); YcRacc_B2(ii)=y(1);
            xRaccR1_B2(ii)=xtan1(1);
            yRaccR1_B2(ii)=ytan1(1);
            xRaccR2_B2(ii)=xtan2(1);
            yRaccR2_B2(ii)=ytan2(1);
        end
        %%
        
        angT2=atan2(yRaccR2_B2(ii)-YcRacc_B2(ii),xRaccR2_B2(ii)-XcRacc_B2(ii));
        angT1=atan2(yRaccR1_B2(ii)-YcRacc_B2(ii),XcRacc_B2(ii)-xRaccR1_B2(ii))+3.14;
        if (angT2<0)
            angT2=pi+angT2;
        end
        
        clear XCer YCer
        XCer=XcRacc_B2(ii)+R_RaccB2(ii)*cos(linspace(angT1,angT2,20));
        YCer=YcRacc_B2(ii)+R_RaccB2(ii)*sin(linspace(angT1,angT2,20));
        %     figure(100);hold on;plot(XCer,YCer,'+m');hold off;
    else
        R_RaccB2(ii)=0;
    end
    
end

%% Cross section drawing
%%
% figure(100);hold on;
% % plot(xo',yo','--r','LineWidth',2); axis([0 rlim 0 rlim]); axis square
% plot(B1k,0,'ob');plot(B2k,0,'ob');
% plot(xpont,ypont,'*c');
% plot(XpBar2,YpBar2,'bs');
% plot(XpBar1,YpBar1,'bs');
% plot(xTraf1,yTraf1,'*m');
% plot(xTraf2,yTraf2,'*m');
% plot(xD1k,yD1k,'ob');
%
% for ii=1:nlay
%    plot([XpBar1(ii),xxD1k(ii)],[YpBar1(ii),yyD1k(ii)],'b','LineWidth',2);
%    plot([XpBar2(ii),xxD2k(ii)],[YpBar2(ii),yyD2k(ii)],'b','LineWidth',2);
%    plot([B1k(ii),XpBar1(ii)],[0,YpBar1(ii)],'b','LineWidth',2);
%    plot([B2k(ii),XpBar2(ii)],[0,YpBar2(ii)],'b','LineWidth',2);
%
% end
% hold off; axis equal;
% keyboard
if (DTrasl~=0)
    a=pi/2-pi/2/p-atan((YpBar1-YpBar2)./(XpBar2-XpBar1));
    Lseg=abs((YpBar1-YpBar2)+1j*(XpBar2-XpBar1));
    hcp=Lseg.*cos(a);
    af=pi/2-pi/2/p-atan((YpBar2(2:end)-YpBar1(1:end-1))./(XpBar1(1:end-1)-XpBar2(2:end)));
    L_f_seg=abs((XpBar1(1:end-1)-XpBar2(2:end))+1j*(YpBar2(2:end)-YpBar1(1:end-1)));
    % hf=[0,L_f_seg.*cos(af)];
    hf=[0,B1k(1:end-1)-B2k(2:end)];
    % hcp=hc;
    for kk=1:nlay
        if kk==1
            B2k(1)=XpBar2(1)+DTrasl;
            [atmp,btmp,ctmp]=retta_per_2pti(B2k(1),0,XpBar2(kk),YpBar2(kk));
            lk(kk)=hcp(kk)/cos(atan(btmp/atmp));
            B1k(1)=B2k(1)-lk(1);
        else
            B2k(kk)=B1k(kk-1)-hf(kk)/cos(atan(btmp/atmp));
            [atmp,btmp,ctmp]=retta_per_2pti(B2k(kk),0,XpBar2(kk),YpBar2(kk));
            lk(kk)=hcp(kk)/cos(atan(btmp/atmp));
            B1k(kk)=B2k(kk)-lk(kk);
        end
        
        %     keyboard
        a1=-atmp/btmp;
        b1=-1;
        c1=(0-a1*B1k(kk));
        m2=tan(pi/2/p);
        a2=m2;
        b2=-1;
        c2=(yTraf1(kk)-m2*xTraf1(kk));
        [x,y]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
        XpBar1(kk)=x;
        YpBar1(kk)=y;
        
    end
    % m1=-a1/b1;
    % q1=-c1/b1;
    % xv=50:01:80;
    % y2=m2*xv+c2;
    % y1=m1*xv+q1;
    % figure(100);plot(xv,y2);hold on;plot(xv,y1); hold off;
end
%%
%%
%%
% Determinazione dei punti caratteristici per il calcolo delle perdite nel
% ferro:
calcolo_posizioni_Pfe;

%% Valutazione ponticelli radiali:
rTemp=abs((x0-xpont)+1j*ypont);
calc_ribs_rad;

if YpBar2(1)<=0
    LowDimBarrier(1)=1;
end

if (LowDimBarrier(1)==1)
    temp.xc=(xxD1k(2:end)+xxD2k(2:end))/2;
    temp.yc=(yyD1k(2:end)+yyD2k(2:end))/2;
    
else
    temp.xc=(xxD1k+xxD2k)/2;
    temp.yc=(yyD1k+yyD2k)/2;
end

% determination of different magnet segment, central point and
% magnetization direction
MagnetFullFill_Seg;

geo.Br = [Br Br];   % doubles Br pieces (half pole + half pole)

temp.xc=xc;
temp.yc=yc;
temp.xmag=xmag;
temp.ymag=ymag;
temp.zmag=zmag;
% Additional division for magnet insertion in flux barrier
temp.XpMag2B1=XpMag2B1;
temp.YpMag2B1=YpMag2B1;
temp.XpMag1B1=XpMag1B1;
temp.YpMag1B1=YpMag1B1;

% end

%%
%% Salvataggio dei dati finali:
temp.B1k=B1k;
temp.B2k=B2k;
temp.Bx0=Bx0;
temp.xpont=xpont;
temp.ypont=ypont;
temp.xxD1k=xxD1k;
temp.yyD1k=yyD1k;
temp.xxD2k=xxD2k;
temp.yyD2k=yyD2k;
temp.XpBar1=XpBar1;
temp.YpBar1=YpBar1;
temp.XpBar2=XpBar2;
temp.YpBar2=YpBar2;

% Coordinate per raccordi barriere di flux:
temp.R_RaccB1=R_RaccB1;
temp.R_RaccB2=R_RaccB2;

temp.XcRacc_B1=XcRacc_B1;
temp.YcRacc_B1=YcRacc_B1;
temp.xRaccR1_B1=xRaccR1_B1;
temp.yRaccR1_B1=yRaccR1_B1;
temp.xRaccR2_B1=xRaccR2_B1;
temp.yRaccR2_B1=yRaccR2_B1;
%
temp.XcRacc_B2=XcRacc_B2;
temp.YcRacc_B2=YcRacc_B2;
temp.xRaccR1_B2=xRaccR1_B2;
temp.yRaccR1_B2=yRaccR1_B2;
temp.xRaccR2_B2=xRaccR2_B2;
temp.yRaccR2_B2=yRaccR2_B2;

temp.xTraf1=xTraf1;
temp.xTraf2=xTraf2;
temp.yTraf1=yTraf1;
temp.yTraf2=yTraf2;

temp.arcLayTraf1=arcLayTraf1;
temp.arcLayTraf2=arcLayTraf2;

temp.XcRibTraf1=XcRibTraf1;
temp.YcRibTraf1=YcRibTraf1;
temp.XcRibTraf2=XcRibTraf2;
temp.YcRibTraf2=YcRibTraf2;

%% Points for radial ribs
temp.XpontRadDx=XpontRadDx;
temp.YpontRadDx=YpontRadDx;
temp.XpontRadSx=XpontRadSx;
temp.YpontRadSx=YpontRadSx;
temp.XpontRadBarDx=XpontRadBarDx;
temp.XpontRadBarSx=XpontRadBarSx;
temp.YpontRadBarDx=YpontRadBarDx;
temp.YpontRadBarSx=YpontRadBarSx;

% add control if the constraction of barrier is too small, so barrier is
% constructed like diamond
temp.LowDimBarrier=LowDimBarrier;
hf=[r,B1k]-[B2k,Ar]; %calcolo dei Delta fi ferro di rotore
geo.hf = hf;
geo.pont = pont;
geo.hc=hc;
geo.xpont=xpont;
geo.ypont=ypont;

% barrier transverse dimension (for permeance evaluation)
temp1_sk = calc_distanza_punti([mean([xxD1k' xxD2k'],2) mean([yyD1k' yyD2k'],2)],[mean([XpBar1' XpBar2'],2) mean([YpBar1' YpBar2'],2)]);
temp2_sk = calc_distanza_punti([mean([B1k' B2k'],2) mean([B1k' B2k'],2)*0],[mean([XpBar1' XpBar2'],2) mean([YpBar1' YpBar2'],2)]);
sk = temp1_sk' + temp2_sk' + 1/0.5 * hc/2;

geo.sk = sk;
geo.pbk = geo.sk ./ geo.hc;
geo.la = sum(geo.hc)/geo.r;
geo.lfe = sum(geo.hf)/geo.r;
geo.ly = (geo.R - (geo.r + geo.g + geo.lt))/geo.r;
geo.B1k=B1k;
geo.B2k=B2k;
geo.hc=hc;
% keyboard
end



