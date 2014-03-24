%% 01/05/2013 Disegno andamento barriere di flusso secondo l'andamento
%% della trasformata di Jukouski w=z^2+(a^2/z)^2 svilupata in coordinate
%% polari...i ribs al traferro sono arrotondati mediante tang alle 2
%% circonferenze generatrici.

function [geo,temp]=nodes_rotor_Fluid(geo)
%%
%% INIZIALIZZAZIONE DATI DI INGRESSO:
%%
% dati0;
l=geo.l;
xr = geo.xr;            % Raggio del rotore al traferro
rlim=xr;
p = geo.p;              % Paia poli
nlay = geo.nlay;        % N� layers
R = geo.r;              % Raggio ext
g = geo.g;              % Traferro
lt = geo.lt;            % Lunghezza denti
pont0 = geo.pont0;      % Ponticelli al traferro
rshaft=geo.Ar;          % raggio albero
Ar=geo.Ar;
nmax=geo.nmax;          % velocit� massima
BandRib=[2.6*geo.pont0,5*geo.pont0,5*geo.pont0];   
rCirRib=1*pont0;
racc_pont=geo.racc_pont;

x0 = xr/cos(pi/2/p);                            % centro cerchi (p generico)
geo.x0 = x0;
Dfe=geo.Dfe;
hc=geo.hc;
alpha=geo.alpha;    %mech angle

alphaRot=geo.alpha*pi/180+pi/2/p;   % angoli barriere di flux al traferro meccanici

%%
%% Determinazione delle costanti per il disegno delle linee mediane in
%%
% corrispondenza degli alphaRot:
C=sin(p*alphaRot).*((((xr-pont0)./rshaft).^(2*p)-1)./((xr-pont0)./rshaft).^p);

Dteta=1;
teta=[0:Dteta*pi/180:pi/p];
r=zeros(length(C),length(teta));
z=zeros(length(C),length(teta));
B0=zeros(length(C),1);

for k=1:length(C)
    
r(k,:)=rshaft*((C(k)+sqrt(C(k)^2+4*(sin(p*teta)).^2))./(2*sin(p*teta))).^(1/p);
z(k,:)=r(k,:).*exp(1j*teta);
Br0(k,:)=rshaft*((C(k)+sqrt(C(k)^2+4))./2).^(1/p);
% coordinata angolare a partire dall'asse y centrale meno il ponticello.
tetaRpont(k)=180/p-(1/p)*asin((C(k)*((xr-geo.pont0)/rshaft)^p)/(((xr-geo.pont0)/rshaft)^(2*p)-1))*180/pi; 
tetaRpont1(k)=180/p-(1/p)*asin((C(k)*((xr-geo.pont0-rCirRib)/rshaft)^p)/(((xr-geo.pont0-rCirRib)/rshaft)^(2*p)-1))*180/pi; 

% % tetaRpont3 sono coordinate che non vengono pi� utilizzate nel disegno del
% % lamierino
% if (k==1)
%     tetaRpont3(k)=180/p-(1/p)*asin((C(k)*((xr-geo.pont0-BandRib1)/rshaft)^p)/(((xr-geo.pont0-BandRib1)/rshaft)^(2*p)-1))*180/pi;
% else
%     tetaRpont3(k)=180/p-(1/p)*asin((C(k)*((xr-geo.pont0-BandRib)/rshaft)^p)/(((xr-geo.pont0-BandRib)/rshaft)^(2*p)-1))*180/pi;
% end
end

x=real(z); y=imag(z);
[xo,yo]=rot_point(x',y',-pi/2/p);

% coordinate in x,y al traferro meno il ponticello
[xpont,ypont]=rot_point((xr-geo.pont0)*cos(tetaRpont*pi/180),(xr-geo.pont0)*sin(tetaRpont*pi/180),-pi/2/p); 
% %pre coordiante in x,y al trefarro dove si trovano D1,2k
% [xpont3,ypont3]=rot_point((xr-geo.pont0-BandRib)*cos(tetaRpont3*pi/180),(xr-geo.pont0-BandRib)*sin(tetaRpont3*pi/180),-pi/2/p); 

%% Fromulazione alternativa per il calcolo della posizione dei centri dei cerchi al traferro

% xp3=xpont-xpont3; yp3=ypont-ypont3;
% tetap3=atan((ypont-ypont3)./(xpont-xpont3));
% D=xp3.*cos(tetap3)+yp3.*sin(tetap3);  E=xp3.^2+yp3.^2-rCirRib^2;
% DeltaRpont=+D-sqrt(D.^2-E)
% xpont1=xpont3+DeltaRpont.*cos(tetap3)
% ypont1=ypont3+DeltaRpont.*sin(tetap3)
%%
% % Old Version of xpont1, ypont1 calc
% % coordinate in x,y al traferro del centro del cerchio di raccordo del Rib
[xpont1,ypont1]=rot_point((xr-geo.pont0-rCirRib)*cos(tetaRpont1*pi/180),(xr-geo.pont0-rCirRib)*sin(tetaRpont1*pi/180),-pi/2/p); 
% rCirRib= sqrt((xpont-xpont1).^2+(ypont1-ypont).^2);
% Valori calcolati ma non significativi per il disegno
% [xpont2,ypont2]=rot_point((xr-geo.pont0-2*rCirRib).*cos(tetaRpont*pi/180),(xr-geo.pont0-2*rCirRib).*sin(tetaRpont*pi/180),-pi/2/p); 

% rotazione di -pi/2/p (rotazione nel 1� e 4� quadrante):
[Bx0,By0]=rot_point(Br0*cos(pi/2/p),Br0*sin(pi/2/p),-pi/2/p);
Bx0=Bx0';
By0=By0';
B1k=Bx0-hc/2+Dfe.*hc/2; B2k=Bx0+hc/2+Dfe.*hc/2;
%%
%% CONTROLLO DI SUCUREZZA PER IL DISEGNO DELLE BARRIERE (NO SOVRAPPOSIZIONI,... RISP. AMPIEZZE MINIME ecc...)
%%
% hc=[1.1992 1.000 3.477];
% Dfe=[-0.2497 -0.75 -0.522];
error_code=[];
Dfe_old=Dfe;
for k=1:nlay  
   if  (hc(k)/2*(1-abs(Dfe(k)))<=pont0);
   disp('#1 Dfe riassegnati');
   Dfe1=1-2*pont0/hc(k);
   error_code=1;
      if (Dfe(k)>0)
       Dfe(k)=Dfe1;
      else
          Dfe(k)=-Dfe1;
      end
   end
end

geo.Dfe=Dfe;
B1k=Bx0-hc/2+Dfe.*hc/2; B2k=Bx0+hc/2+Dfe.*hc/2;

if nlay~=1
    
for k=1:nlay-1
    
    hfemin=2*pont0;
    hcmin=1;
    hc_old=hc;
    if ((xr-B2k(k))<1) % questa condizione varrebbe per la 1�barriera
        B2k(k)=xr-1;
        disp('#2 vincolo 1 layer esce dal rotore');
        error_code=[error_code,2];
    end
    if (B1k(end)<geo.Ar+hfemin)    % questa condizione vale per l'ultima barriera di flux
        B1k(end)=geo.Ar+1.0;
        disp('#3 vincolo n� layer interseca albero')   
        error_code=[error_code,3];
    end                

    if (B2k(k+1)>=B1k(k))   % questa condizione vale invece per tutte le barriere
        Dq=B2k(k+1)-B1k(k);
        B2p=B2k(k+1)-(1/2)*(Dq+hfemin);
        B1p=B1k(k)+(1/2)*(Dq+hfemin);
        disp('#4 vincolo 1-n� layer sovrapposizione'); 
        error_code=[error_code,4];
       
        if ((Bx0(k)<B1p)||(Bx0(k+1)>B2p)|| ((B2p-Bx0(k+1))<pont0/2) || ((Bx0(k)-B1p)<pont0/2))  % condizione vale nel caso in cui muovendosi con #2 non c'� pi� spazio per l'aria
           B1p=Bx0(k)-(Bx0(k)-Bx0(k+1))/3;
           B2p=Bx0(k+1)+(Bx0(k)-Bx0(k+1))/3;
           disp('#5 vincolo 1-n� intersezione arie, spessore lato barriera<pont0/2 --> equa ripartizione aria ferro');
           error_code=[error_code,5];
        end
        B2k(k+1)=B2p;
        B1k(k)=B1p;
        
    elseif((B1k(k)-B2k(k+1))<hfemin)    %questa condizione vale quando non si � riusciti ad assicurare un ferro minimo tra le barriere di flux (non ho trovato nulla di meglio per adesso, pensarci su!!!!)
        dB12=Bx0(k)-Bx0(k+1);
        Dhc12=dB12-hfemin;
        B1p=Bx0(k)-Dhc12/2;
        B2p=Bx0(k+1)+Dhc12/2;
        disp('#6 vincolo 1-n� tra 2 barriere successive no abbastanza ferro --> equa ripartizione aria ferro');        
        error_code=[error_code,6];
        
        if ((Bx0(k)-B1p)<pont0 || (B2p-Bx0(k+1))<pont0)
        B1p=Bx0(k)-(Bx0(k)-Bx0(k+1))/3;
        B2p=Bx0(k+1)+(Bx0(k)-Bx0(k+1))/3;
        disp('#7 vincolo 1-n� (Bx0(k)-B1p)<pont0 (B2p-Bx0(k+1))<pont0 --> equa ripartizione aria ferro');
        error_code=[error_code,7];        
        end
        
        B2k(k+1)=B2p;
        B1k(k)=B1p;

    end
%
%     if ((Bx0(k)-B1k(k))<(pont0/2))
%         B1k(k)=Bx0(k)+pont0/2;
%                 
%     end
    
end

else
    hfemin=2*pont0;
    hcmin=1;
    hc_old=hc;
    if ((xr-B2k(k))<1) % questa condizione varrebbe per la 1�barriera
        B2k(k)=xr-1;
        disp('#2 vincolo 1 layer esce dal rotore');
        error_code=[error_code,2];
    end

end
error_code=[error_code,zeros(1,3*nlay-length(error_code))];

hc=B2k-B1k; geo.hc=hc;

% Instruction for flux barrier translaction in function of range Dfe=[-1:1]:
for k=1:length(Dfe)
if (Dfe(k)==1)||(Dfe(k)==-1)
    Bx0(k)=(B1k(k)+B2k(k))/2;
    [xBx0New(k),yBy0New(k)]=rot_point(Bx0(k),0,pi/2/p);
    Br0New(k)=abs(xBx0New(k)+1j*yBy0New(k));
    CNew(k)=(((Br0New(k)./rshaft).^(2*p)-1)./(Br0New(k)./rshaft).^p);
    tetaRpontNew(k)=180/p-(1/p)*asin((CNew(k)*((xr-geo.pont0)/rshaft)^p)/(((xr-geo.pont0)/rshaft)^(2*p)-1))*180/pi; 
    % Riassegno xpont ed ypont
    [xpontNew,ypontNew]=rot_point((xr-geo.pont0)*cos(tetaRpontNew(k)*pi/180),(xr-geo.pont0)*sin(tetaRpontNew(k)*pi/180),-pi/2/p); 
    xpont(k)=xpontNew;
    ypont(k)=ypontNew;
end
end

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%
%% SECURITY CONTROL END
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%

%% ROTOR LAYER DRAWING:
% si ruotano i valori di Bk nel 1� quadrante ( poich� la validit� delle formule per la forma delle barriere � limitata al 1� quadrante):

[xB1k,yB1k]=rot_point(B1k,0,pi/2/p);
[xB2k,yB2k]=rot_point(B2k,0,pi/2/p);

rB1k=abs(xB1k+1j*yB1k);
rB2k=abs(xB2k+1j*yB2k);
CB1k=(((rB1k./rshaft).^(2*p)-1)./(rB1k./rshaft).^p);
CB2k=(((rB2k./rshaft).^(2*p)-1)./(rB2k./rshaft).^p);

% Ciclo per il disegno dei lati delle barriere di flux:
for k=1:length(CB1k)
    
RB1k(k,:)=rshaft*((CB1k(k)+sqrt(CB1k(k)^2+4*(sin(p*teta)).^2))./(2*sin(p*teta))).^(1/p);
zB1k(k,:)=RB1k(k,:).*exp(1j*teta);

RB2k(k,:)=rshaft*((CB2k(k)+sqrt(CB2k(k)^2+4*(sin(p*teta)).^2))./(2*sin(p*teta))).^(1/p);
zB2k(k,:)=RB2k(k,:).*exp(1j*teta);

tetaIBk1(k)=180/p-(1/p)*asin((CB1k(k)*((xr-geo.pont0)/rshaft)^p)/(((xr-geo.pont0)/rshaft)^(2*p)-1))*180/pi;
tetaIBk2(k)=180/p-(1/p)*asin((CB2k(k)*((xr-geo.pont0)/rshaft)^p)/(((xr-geo.pont0)/rshaft)^(2*p)-1))*180/pi;
tetaIBK2=45;
tetaTraf01(k)=180/p-(1/p)*asin((CB1k(k)*(xr/rshaft)^p)/((xr/rshaft)^(2*p)-1))*180/pi;
tetaTraf02(k)=180/p-(1/p)*asin((CB2k(k)*(xr/rshaft)^p)/((xr/rshaft)^(2*p)-1))*180/pi;

end

% rotazione dei valori nel 1�e 4� quadrante pronti per il disegno in matlab
% del lamierino (no FEMM)
[xxB1k,yyB1k]=rot_point(real(zB1k),imag(zB1k),-pi/2/p);
[xxB2k,yyB2k]=rot_point(real(zB2k),imag(zB2k),-pi/2/p);

% Determinazione xTraf e yTraf punti al traferro:
[xTraf01,yTraf01]=rot_point(xr*cos(tetaTraf01*pi/180),xr*sin(tetaTraf01*pi/180),-pi/2/p);
[xTraf02,yTraf02]=rot_point(xr*cos(tetaTraf02*pi/180),xr*sin(tetaTraf02*pi/180),-pi/2/p);
% Determinazione xTraf e yTraf punti al traferro - pont0:
[xTraf1,yTraf1]=rot_point((xr-geo.pont0)*cos(tetaIBk1*pi/180),(xr-geo.pont0)*sin(tetaIBk1*pi/180),-pi/2/p);
[xTraf2,yTraf2]=rot_point((xr-geo.pont0)*cos(tetaIBk2*pi/180),(xr-geo.pont0)*sin(tetaIBk2*pi/180),-pi/2/p);

% Riassegnazione ed interpolazione dei lati di barriera su 50 punti solo
% per la parte positiva in xy di geo mot
[x,y]=interp_flux_barrier(xxB1k,yyB1k,2*yTraf01);
clear xxB1k yyB1k;
xxB1k=x; yyB1k=y;
clear x y
[x,y]=interp_flux_barrier(xxB2k,yyB2k,2*yTraf02);
clear xxB2k yyB2k;
xxB2k=x; yyB2k=y;

%%
%% DISEGNO DEI PUNTI AL TRAFERRO MEDIANTE TANGENTE TRA CIRCONFERENZE:
%%

for k=1:length(CB2k)
    [xxB2k_med,yyB2k_med]=valore_medio_di_barriera(xxB2k(k,:),yyB2k(k,:),yTraf2(k));
    [x02,y02,r02]=circonferenza_per_3_pti(B2k(k),0,xTraf2(k),yTraf2(k),xxB2k_med,yyB2k_med);
    r=xr-pont0;
    [xD2k,yD2k,xc,yc,rc]=cir_tg_2cir(xpont(k),ypont(k),r,x02,y02,r02);
    angT1=atan2(yD2k-yc,xD2k-xc);
    angT2=atan2(ypont(k)-yc,xpont(k)-xc);
%     if (angT1<angT2)
%         angT2=-2*pi+angT2;
%     end
    XcRibTraf2(k)=xc;
    YcRibTraf2(k)=yc;
    arcLayTraf2(k)=(angT2-angT1)*180/pi;
    xxD2k(k)=xD2k;
    yyD2k(k)=yD2k;
        clear XCer YCer
    XCer=xc+rc*cos([angT1:0.01:angT2]);
    YCer=yc+rc*sin([angT1:0.01:angT2]);
%     figure(100);plot(XCer,YCer,'k');hold on;
end


for k=1:length(CB1k)
    [xxB1k_med,yyB1k_med]=valore_medio_di_barriera(xxB1k(k,:),yyB1k(k,:),yTraf1(k));
    [x02,y02,r02,angleA,angleB]=circonferenza_per_3_pti(B1k(k),0,xTraf1(k),yTraf1(k),xxB1k_med,yyB1k_med);
    if (angleA<angleB)
        angleB=-2*pi+angleB;
    end
    XCerchio=x02+r02*cos([angleB:0.01:angleA]);
    YCerchio=y02+r02*sin([angleB:0.01:angleA]);
%     figure(100);hold on;plot(XCerchio,YCerchio,'k');hold off;
    
    r=xr-pont0;
    [xD1k,yD1k,xc,yc,rc]=cir_tg_2cir(xpont(k),ypont(k),r,x02,y02,r02);
    angT2=atan2(yD1k-yc,xD1k-xc);
    angT1=atan2(ypont(k)-yc,xpont(k)-xc);
%     figure(100);hold on;plot(xD1k,yD1k,'rs');plot(xc,yc,'gs'); hold off;
%     figure(100);hold on;plot([0,xpont(k)]',[0,ypont(k)]','-b'); hold off;

    if (angT2<0)
        angT2=pi+angT2;
    end
    arcLayTraf1(k)=(angT2-angT1)*180/pi;
    XcRibTraf1(k)=xc;
    YcRibTraf1(k)=yc;
    xxD1k(k)=xD1k;
    yyD1k(k)=yD1k;
    clear XCer YCer
    XCer=xc+rc*cos(linspace(angT1,angT2,20));
    YCer=yc+rc*sin(linspace(angT1,angT2,20));
%     figure(100);hold on;plot(XCer,YCer,'k');hold off;
    

end

clear xxB1k_mean yyB1k_mean xxB2k_mean yyB2k_mean;

% figure(100);hold on; 
% plot(xo,yo,'--r','LineWidth',2); axis([0 rlim 0 rlim]); axis square
% % plot(xr*cos(geo.alpha*pi/180),xr*sin(geo.alpha*pi/180),'or');
% plot(xpont,ypont,'or');
% plot(xTraf1,yTraf1,'ob',xTraf2,yTraf2,'ob');
% plot(xxD1k,yyD1k,'bs');
% plot(xxD2k,yyD2k,'bs');
% xx=(0:0.1:rlim);
% plot(xx,sqrt(rlim^2-xx.^2),'k','LineWidth',2);
% plot([0:0.1:rshaft],sqrt(rshaft^2-[0:0.1:rshaft].^2),'k','LineWidth',2);
% plot(xxB1k',yyB1k','b','LineWidth',2);plot(xxB2k',yyB2k','b','LineWidth',2);
% plot(B1k,0,'ob');plot(B2k,0,'ob');
% % plot(xRib1,yRib1,'ms',xRib2,yRib2,'ms');
% % plot(xpont1,ypont1,'om');
% hold off;grid on;
%%
%% SET DI ISTRUZIONI INTERPOLAZIONE BARRIERE DI FLUX:
%%
[xxB1k_mean,yyB1k_mean]=valore_medio_di_barriera(xxB1k,yyB1k,yTraf1);
[xxB2k_mean,yyB2k_mean]=valore_medio_di_barriera(xxB2k,yyB2k,yTraf2);
% figure(100);hold on; plot(xxB1k_mean,yyB1k_mean,'og',xxB2k_mean,yyB2k_mean,'om');hold off;
%%

for i=1:nlay
    if (i==1)
        % if (B2k(i)>=(xr-pont0-BandRib(i))||yyB2k(i,2)>=yyD2k(i)||yyD2k(1)==0||isnan(yyD2k(1)))
        if (yyB2k(i,2)>=yyD2k(i)||yyD2k(1)<=0||isnan(yyD2k(1)))
            case1layer=1;
            tetaBuccia2k=180/p-(1/p)*asin((CB2k(i)*((xr-geo.pont0-BandRib(i))/rshaft)^p)/(((xr-geo.pont0-BandRib(i))/rshaft)^(2*p)-1))*180/pi;
            % rotazione nel 1� e 4� quadrante:
            [xBuccia2k,yBuccia2k]=rot_point((xr-geo.pont0-BandRib(i))*cos(tetaBuccia2k*pi/180),(xr-geo.pont0-BandRib(i))*sin(tetaBuccia2k*pi/180),-pi/2/p);
            xxD2k(1)=xBuccia2k;
            yyD2k(1)=yBuccia2k;
        else
            case1layer=0;
            [xxB2k_mean2,yyB2k_mean2]=valore_medio_di_barriera(xxB2k(i,:),yyB2k(i,:),yyD2k(i));
            [xc,yc,r,angleA,angleB]=circonferenza_per_3_pti(B2k(i),0,xxD2k(i),yyD2k(i),xxB2k_mean2,yyB2k_mean2);
            if (angleA<angleB)
                angleB=-2*pi+angleB;
            end
            XCerchio=xc+r*cos(linspace(angleB,angleA,20));
            YCerchio=yc+r*sin(linspace(angleB,angleA,20));
            XcBar2(i)=xc;
            YcBar2(i)=yc;
            
%             figure(100);hold on;plot(XCerchio,YCerchio,'k');hold off;
            arcLayer2(i)=(angleA-angleB)*180/pi;
            if abs(arcLayer2(i))>90
                arcLayer2(i)=360-abs(arcLayer2(i));
            end
            
        end
    else
        [xc,yc,r,angleA,angleB]=circonferenza_per_3_pti(B2k(i),0,xxD2k(i),yyD2k(i),xxB2k_mean(i),yyB2k_mean(i));
        
        if (angleA<angleB)
            angleB=-2*pi+angleB;
        end
        XCerchio=xc+r*cos(linspace(angleB,angleA,20));
        YCerchio=yc+r*sin(linspace(angleB,angleA,20));
        XcBar2(i)=xc;
        YcBar2(i)=yc;
%         figure(100);hold on;plot(XCerchio,YCerchio,'k');hold off;
        arcLayer2(i)=(angleA-angleB)*180/pi;
        if abs(arcLayer2(i))>90
            arcLayer2(i)=360-abs(arcLayer2(i));
        end
        
    end
    
end

if nlay~=1
    
for i=1:nlay-1
        [xc,yc,r,angleA,angleB]=circonferenza_per_3_pti(B1k(i),0,xxD1k(i),yyD1k(i),xxB1k_mean(i),yyB1k_mean(i));
        
        if (angleA<angleB)
            angleB=-2*pi+angleB;
        end
        
        XCerchio=xc+r*cos(linspace(angleB,angleA,20));
        YCerchio=yc+r*sin(linspace(angleB,angleA,20));
        
        XcBar1(i)=xc;
        YcBar1(i)=yc;

%     figure(100);hold on;plot(XCerchio,YCerchio,'k');hold off;
    arcLayer1(i)=(angleA-angleB)*180/pi;
end
%%
%% SET ISTRUZIONI INTERPOLAZIONI ULTIMA BARRIERA DI FLUX:
%%
% index=find(yyB1k(nlay,:)>=(yyB1k_mean(end)+yTraf1(end))/2);
% yyB1k_mean2=yyB1k(nlay,index(1));    
% xxB1k_mean2=xxB1k(nlay,index(1));

[xxB1k_mean2,yyB1k_mean2]=valore_medio_di_barriera(xxB1k(nlay,:),yyB1k(nlay,:),yyB1k_mean(end));

% figure(100);hold on; plot(xxB1k_mean2,yyB1k_mean2,'bs');hold off;
clear index;

% index=find(yyB1k(nlay,:)>=yyB1k_mean2 & yyB1k(nlay,:)<=(yyB1k_mean(end)));
% yyB1kAux3=yyB1k(nlay,index);    
% xxB1kAux3=xxB1k(nlay,index);
% [xxB1k_mean3,yyB1k_mean3]=valore_medio_di_barriera(xxB1kAux3,yyB1kAux3,yyB1k_mean(end));

%% 2013/10/04 MG determinazione del segmento intermedio

index=find(yyB1k(nlay,:)>=(yyB1k_mean(end)) & yyB1k(nlay,:)<=0.9*yyD1k(end));

if length(index)==1
    index=find(yyB1k(nlay,:)>=(yyB1k_mean(end)) & yyB1k(nlay,:)<=1.1*yyD1k(end));
end

yyB1kAux3=yyB1k(nlay,index);
xxB1kAux3=xxB1k(nlay,index);
[xxB1k_mean3,yyB1k_mean3]=valore_medio_di_barriera(xxB1kAux3,yyB1kAux3,0.9*yyD1k(end));

% figure(100);hold on; plot(xxB1k_mean3,yyB1k_mean3,'k');hold off;
clear index;

% [xc,yc,r,angleA,angleB]=circonferenza_per_3_pti(xxB1k_mean2,yyB1k_mean2,xxB1k_mean(end),yyB1k_mean(end),xxB1k_mean3,yyB1k_mean3);
[xc,yc,r,angleA,angleB]=circonferenza_per_3_pti(xxB1k_mean2,yyB1k_mean2,xxB1k_mean3,yyB1k_mean3,xxB1k_mean(end),yyB1k_mean(end));
% centro raggio ultima barriera di flux:
    XcBarLast_mean=xc;
    YcBarLast_mean=yc;
    
    if (angleA<angleB)
        angleB=-2*pi+angleB;
    end
    XCerchio=xc+r*cos([ angleB:0.01:angleA]);
    YCerchio=yc+r*sin([ angleB:0.01:angleA]);
    
%     figure(100);hold on;plot(XCerchio,YCerchio,'k');hold off;
    arcLayer3=(angleA-angleB)*180/pi;
    [xD1k,yD1k,xc,yc,rc]=tg_cir(xxB1k_mean3,yyB1k_mean3,xTraf1(end),yTraf1(end),xpont(end),ypont(end));
    XcRibTraf1(end)=xc;
    YcRibTraf1(end)=yc;
    angT2=atan2(yD1k-yc,xD1k-xc);
    angT1=atan2(ypont(k)-yc,xpont(k)-xc);
%     figure(100);hold on;plot(xD1k,yD1k,'k');plot(xc,yc,'k'); hold off;
%     figure(100);hold on;plot([0,xpont(k)]',[0,ypont(k)]','k'); hold off;

    if (angT2<0)
        angT2=pi+angT2;
    end
    arcLayTraf1(end)=(angT2-angT1)*180/pi;
    xxD1k(end)=xD1k;
    yyD1k(end)=yD1k;
    clear XCer YCer
    XCer=xc+rc*cos(linspace(angT1,angT2,20));
    YCer=yc+rc*sin(linspace(angT1,angT2,20));
%     figure(100);hold on;plot(XCer,YCer,'k');hold off;

else
 
        [xc,yc,r,angleA,angleB]=circonferenza_per_3_pti(B1k(1),0,xxD1k(1),yyD1k(1),xxB1k_mean(1),yyB1k_mean(1));
        
        if (angleA<angleB)
            angleB=-2*pi+angleB;
        end
        
        XCerchio=xc+r*cos(linspace(angleB,angleA,20));
        YCerchio=yc+r*sin(linspace(angleB,angleA,20));
        
        XcBar1(1)=xc;
        YcBar1(1)=yc;

%     figure(100);hold on;plot(XCerchio,YCerchio,'k');hold off;
    arcLayer1(1)=(angleA-angleB)*180/pi;


xxB1k_mean2=NaN; 
yyB1k_mean2=NaN;
xxB1k_mean3=NaN; 
yyB1k_mean3=NaN;
XcBarLast_mean=NaN;
YcBarLast_mean=NaN;
arcLayer3=NaN;

end

%%
%%
%% Calcolo dei volumi di ferro:
%%
calc_Iron_area_fluid;
%%   
% Determinazione dei punti caratteristici per il calcolo delle perdite nel
% ferro:
calcolo_posizioni_Pfe;
%% Salvataggio dei dati finali:
%%
hf=[xr,B1k]-[B2k,Ar]; % rotor Dfe calculation in absolute term.
%% variabili da portare fuori 
geo.hf = hf;

temp.xpont1=xpont1;
temp.ypont1=ypont1;

temp.Bx0=Bx0;
temp.By0=By0;
temp.B1k=B1k;
temp.B2k=B2k;

temp.xpont=xpont;
temp.ypont=ypont;

% temp.RcirRib=rCirRib;
% geo.xpont2=xpont2;
% geo.ypont2=ypont2;

temp.xTraf1=xTraf1;
temp.xTraf2=xTraf2;
temp.yTraf1=yTraf1;
temp.yTraf2=yTraf2;
temp.arcLayer1=arcLayer1;
temp.arcLayer2=arcLayer2;
temp.arcLayTraf1=arcLayTraf1;
temp.arcLayTraf2=arcLayTraf2;
temp.XcRibTraf1=XcRibTraf1;
temp.YcRibTraf1=YcRibTraf1;
temp.XcRibTraf2=XcRibTraf2;
temp.YcRibTraf2=YcRibTraf2;
%% Punti per i ribs radiali
temp.XpontRadSx=XpontRadSx;
temp.YpontRadSx=YpontRadSx;
temp.XpontRadDx=XpontRadDx;
temp.YpontRadDx=YpontRadDx;
temp.XpontRadBarDx=XpontRadBarDx;
temp.YpontRadBarDx=YpontRadBarDx;
temp.XpontRadBarSx=XpontRadBarSx;
temp.YpontRadBarSx=YpontRadBarSx;

temp.xxD1k=xxD1k;
temp.yyD1k=yyD1k;
temp.xxD2k=xxD2k;
temp.yyD2k=yyD2k;

temp.XcBar1=XcBar1;
temp.YcBar1=YcBar1;
temp.XcBar2=XcBar2;
temp.YcBar2=YcBar2;
temp.XcBarLast_mean=XcBarLast_mean;
temp.YcBarLast_mean=YcBarLast_mean;

temp.xxB1k_mean=xxB1k_mean;
temp.yyB1k_mean=yyB1k_mean;
temp.xxB2k_mean=xxB2k_mean;
temp.yyB2k_mean=yyB2k_mean;
temp.xxB1k_mean2=xxB1k_mean2;
temp.yyB1k_mean2=yyB1k_mean2;
temp.xxB1k_mean3=xxB1k_mean3;
temp.yyB1k_mean3=yyB1k_mean3;
temp.arcLayer3=arcLayer3;
temp.error_code=error_code;
temp.r_fe=rfe;
temp.x_fe=xFe;
temp.y_fe=yFe;
geo.pont=pont;
% keyboard

end

