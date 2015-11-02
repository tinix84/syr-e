%% DISEGNO DEL ROTORE (GRUPPO 0)
%% modifica 15 novembre - gp
% - corretto angolo di magnetizzazione errato
% - aggiunto albero per evitare saturazione spider (irrealistica)
%% modifica 08 gennaio 2010 - gp
% - porto fuori geo.r_all (servono per Pfe)
% - calcolo e porto fuori geo.r_fe,geo.hf (servono per Pfe)

function [geo] = disegna_rotore_semicerchi(geo,mat)
 
global xy_ferro

group0 = 0;                      % di default il rotore è gruppo 0, alla fine dello script gli assegno il gruppo 2 (mod. 28 04 10)

% NOMENCLATURA: PROPRIETA' GEOMETRICHE E MATERIALI

xr = geo.xr;                    % Raggio del rotore al traferro
x0 = geo.x0;                    % Centro fittizio
Ar = geo.Ar;                    % Raggio albero
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

magnet = mat.magnet;            % Materiale magneti: 'NdFeB 37 MGOe'
steel = mat.steel;              % Materiale ferromagnetico: 'M-27 Steel'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ampiezza di un arco lungo pont0
% ang_pont0 = pont0 / xr * 180/pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CONTORNO DEL ROTORE: INDIVIDUAZIONE DELLE COORDINATE PRINCIPALI; DISEGNO
%% DELL'ALBERO E DELLE LINEE DI INZIO E FINE POLO.
% DEFINIZIONE DELLE COORDINATE CHE INDIVIDUANO IL CONTORNO PRINCIPALE DEL ROTORE VERSO IL TRAFERRO

x1 = xr;        y1 = 0;                             % Punto1
[x2,y2] = rot_point(x1,y1,180/p*pi/180);            % Punto2
x3 = xr-pont0;  y3 = 0;                             % Punto3
[x4,y4] = rot_point(x3,y3,180/p*pi/180);            % Punto4

% LINEA INIZIO POLO (CONDIZIONE AL CONTORNO)
mi_drawline(0,0,x1,y1);                             % mi_drawline(x1,y1,x2,y2): Aggiunge i nodi (x1,y1) e (x2,y2), e disegna una linea tra questi due
% punti.
mi_selectsegment(mean([0 x1]),mean([0 y1]));        % Con mi_selectsegment(x,y) si seleziona la linea più vicina a (x,y)
mi_setsegmentprop('APr1', 0, 1, 0, group0);          % Il comando mi_setsegmentprop(’propname’, elementsize, automesh, hide, group0) impone che i
                                                    % segmenti selezionati abbiano le seguenti proprietà: 1) la condizione al contorno 'propname';
                                                    % 2)e 3)la dimensione dei loro elementi locali deve essere inferiore a elementsize (automesh=0)
                                                    % oppure deve essere scelta automaticamente (automesh=1). 4) se hide=0 i segmenti non devono
                                                    % nascosti nel post-processor; 5) infine il numero corrispondente a group0 indica a quale gruppo
                                                    % devono appartenere i segmenti.

% LINEA FINE POLO (CONDIZIONE AL CONTORNO)

mi_drawline(0,0,x2,y2);
mi_selectsegment(mean([0 x2]),mean([0 y2]));

mi_setsegmentprop('APr1', 0, 1, 0, group0);
mi_clearselected;
%keyboard
% ALBERO (CONDIZIONE AL CONTORNO)

%% modifica 15 novembre: disegno l'albero e assegno Anti Periodicità (APr0)
mi_addboundprop('APr0', 0, 0, 0, 0, 0, 0, 0, 0, 5);

x1 = Ar; y1 = 0;                                    % Punto1
[x2,y2] = rot_point(x1,y1,180/p*pi/180);            % Punto2
mi_drawarc(x1,y1,x2,y2,180/p,10);                   % 
mi_selectsegment(mean([0 x1]),mean([0 y1]));
mi_selectsegment(mean([0 x2]),mean([0 y2]));

mi_setsegmentprop('APr0', 0, 1, 0, group0);
mi_clearselected;

mi_addblocklabel(1,1);
mi_selectlabel(1,1);
mi_setblockprop(steel, 1, 0, 'None', 0, group0, 1);

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

%% CALCOLO PONTICELLI IN BASE ALLA CENTRIFUGA
% 1. volume guide di flusso
% 2. somma dei volumi e calcolo delle masse sospese su ciascun pont 
% 3. calcolo dello spessore pont che corrisponde alla sollecitazione
% centrifuga voluta (sigma_max = 180 N/mm2)

% Volume singolo settore di ferro
% (i settori di ferro in cui si divide il rotore sono separati dalle circonferenze di raggio r centrate in (x0,0))
%dVol_Fe = [];	
%rG = [];                        % Distanza da (0,0) dei baricentri dei volumi appesi ai ponticelli

%for jj = 1:2:length(r_all)
%    kk = ceil(jj/2);            % indice ponticello
    
%    if kk == 1

 %       dVol_Fe(1) = l * (r(1)^2*(beta(1)*pi/180-0.5*sin(2*beta(1)*pi/180)) + xr^2*(alpha(1)*pi/180-0.5*sin(2*alpha(1)*pi/180)));

  %  else

   %     dVol_Fe(kk) = l * (r(kk)^2-r(kk-1)^2)*0.5*(beta(kk)*pi/180+beta(kk-1)*pi/180);

    %end

%end

%Vol_Fe = integra_fx(1:length(dVol_Fe),dVol_Fe);      % volumi appesi ai ponticelli
%M_Fe = Vol_Fe * 1e-9 * 7800 ;                        % massa ferro appeso ai ponticelli

%rG =  0.5 * (xr + (x0 - r));

%F_centrifuga = M_Fe .* rG/1000 *  (nmax * pi/30)^2;
%sigma_max = 180;                                     % N/mm2 - snervamento lamierino

%pont = F_centrifuga/(sigma_max * l);                 % mm

%pont(pont < pont0) = 0;                              % NOTA BENE: Elimino i ponticelli troppo sottili
 %pont(pont < pont0) = pont0;                              % NOTA BENE: non permetto ponticelli troppo sottili (riga già commentata in versione originale)

%% DISEGNO DELLE BARRIERE DI FLUSSO (ASSEGNO MAGNETI E ARIA).
%keyboard
mi_selectgroup(group0); mi_moverotate(0,0,-90/p);                                    % Ruoto il contorno che ho già tracciato (linee di inizio e fine
% polo + albero)in maniera che il disegno sia simmetrico rispetto
% all'asse x.
X4 = zeros(1,nlay); Y4 = X4;
X42=[];Y42=X42; X32=[];Y32=X32;

 jj = 1:2:length(r_all);                                                          
    [xc,yc] = calc_intersezione_cerchi(xr-pont0-hc(ceil(jj/2))/2, r(ceil(jj/2)), x0);
   
[alphacc,rcc] = cart2pol(xc,yc);
alphaccr=alphacc+pi/2/p;
%keyboard

ff=2;
i=1;
for jj=1:1:nlay        %sostituire 3 con nlay è più generale!!
    if (jj>1)
x3=hc(jj)/2+rcc(jj)*cos(alphacc(jj)+pi/2/p);
x4=-hc(jj)/2+rcc(jj)*cos(alphacc(jj)+pi/2/p);
y3=rcc(jj)*sin(alphacc(jj)+pi/2/p); 
y4=y3;
x3n=y3; x4n=y4; y3n=x3; y4n=x4;
%rotazione di -pi/2p e aggiorno le coordinate:
[alphar,rr] = cart2pol(x3,y3);
x3=rr*cos(alphar-pi/2/p);
y3=rr*sin(alphar-pi/2/p);
[alphar,rr] = cart2pol(x4,y4);
x4=rr*cos(alphar-pi/2/p);
y4=rr*sin(alphar-pi/2/p);

[alphar,rr] = cart2pol(x3n,y3n);
x3n=rr*cos(alphar-pi/2/p);
y3n=rr*sin(alphar-pi/2/p);

[alphar,rr] = cart2pol(x4n,y4n);
x4n=rr*cos(alphar-pi/2/p);
y4n=rr*sin(alphar-pi/2/p);

X4(i)=x4;
Y4(i)=y4;
X3(i)=x3;
Y3(i)=y3;
%keyboard

hfe2=r_all(ff+1)-r_all(ff);
ff=jj+2;
x32=x2-hfe2;
a=x3-x32;
y32=y3-a*tan(pi/4);
x42=x32-hc(jj);
d=x4-x42;
y42=y4-d*tan(pi/4);
x2=x42;
Xbar(i)=x42;
X42(i)=x42;
Y42(i)=y42;
X32(i)=x32;
Y32(i)=y32;

if x3<x32 | y3<y32
x3=x32;
y3=y32;
y3n=-y3;
x3n=x3;
else
end

mi_addnode(x3,y3);
mi_addnode(x4,y4);
mi_addnode(x3n,y3n);
mi_addnode(x4n,y4n);
mi_drawarc(x3,y3,x4,y4,180,10);
mi_drawarc(x4n,y4n,x3n,y3n,180,10);
mi_addnode(x32,y32);
mi_addnode(x42,y42);
mi_addnode(x32,-y32);
mi_addnode(x42,-y42);
mi_drawline(x3,y3,x32,y32);
mi_drawline(x4,y4,x42,y42);
mi_drawline(x32,y32,x32,-y32);
mi_drawline(x42,y42,x42,-y42);
mi_drawline(x32,y32,x3,y3);
mi_drawline(x42,y42,x4,y4);
mi_drawline(x32,-y32,x3n,y3n);
mi_drawline(x42,-y42,x4n,y4n);
        x10 =x32-racc_pont;                           % Coordinata x Nodo 10 ((Ponticello radiale_vertice in alto a dx)=(Punto
        % in basso a dx del raccordo))
        x11 = x42+racc_pont;                         % Coordinata x Nodo 11 ((Ponticello radiale_vertice in alto a sx)=(Punto
        % in basso a sx del raccordo))
    else 
        x1=xc(jj)+hc(jj)/2; x2=xc(jj)-hc(jj)/2;
        y1=yc(jj); y2=yc(jj); y1n=-y1; y2n=-y2;
        mi_addnode(x1,y1);mi_addnode(x2,y2);mi_addnode(x1,y1n);mi_addnode(x2,y2n);
        mi_drawarc(x1,y1,x2,y2,180,10); mi_drawarc(x2,y2n,x1,y1n,180,10);
        mi_drawline(x1,y1,x1,y1n);
        mi_drawline(x2,y2,x2,y2n);
        x10 =x1-racc_pont;                           % Coordinata x Nodo 10 ((Ponticello radiale_vertice in alto a dx)=(Punto
        % in basso a dx del raccordo))
        x11 = x2+racc_pont;                         % Coordinata x Nodo 11 ((Ponticello radiale_vertice in alto a sx)=(Punto
        % in basso a sx del raccordo))
        Xbar(i)=x2;
        X4(i)=x2;
        X42(i)=x2;
        Y4(i)=y2;
        Y42(i)=y2;
        X32(i)=x1;
        Y32(i)=y1;

    end

%% (Matteo 16/11/2011 calcolo sforzi centrifughi):
%% CALCOLO PONTICELLI IN BASE ALLA CENTRIFUGA
% 1. volume guide di flusso
% 2. somma dei volumi e calcolo delle masse sospese su ciascun pont 
% 3. calcolo dello spessore pont che corrisponde alla sollecitazione
% centrifuga voluta (sigma_max = 180 N/mm2)

% Volume singolo settore di ferro
% (i settori di ferro in cui si divide il rotore sono separati dalle circonferenze di raggio r centrate in (x0,0))
dVol_Fe = [];	
rG = [];   % Distanza da (0,0) dei baricentri dei volumi appesi ai ponticelli
ii=jj-1;
kk=jj;
pont=[];
    if kk == 1

        dVol_Fe(1) = l*xr^2*(alpha(1)*pi/180-0.5*sin(2*alpha(1)*pi/180));

    else
        area_guide_ferro_rotore;
        dVol_Fe(kk) = l * Atot;

    end

Vol_Fe = integra_fx(1:length(dVol_Fe),dVol_Fe);      % volumi appesi ai ponticelli
Vol_Fe=Vol_Fe(1,kk);
M_Fe = Vol_Fe * 1e-9 * 7800 ;                        % massa ferro appeso ai ponticelli

rG=0.5*(xr+(X32(i)+X42(i))/2);

F_centrifuga = M_Fe .* rG/1000 *  (nmax * pi/30)^2;
sigma_max = 180;                                     % N/mm2 - snervamento lamierino

pont(jj) = F_centrifuga/(sigma_max * l);                % mm

if pont(jj)<pont0
    pont(jj)=0;                             % NOTA BENE: Elimino i ponticelli troppo sottili
else
end
                            
%%

        if (x11 < x10) && (pont(ceil(jj))>0)                   % Se i due punti non si incrociano (raccordo non troppo grande rispetto alla
            % larghezza della barriera), e se lo spessore estratto dall'algoritmo per il
            % ponticello non è troppo piccolo, allora procedo al disegno del ponticello
%keyboard
            lpont = x10 - x11;                                      % lpont=lunghezza ponticello=larghezza barriera-2*racc_pont
            hpont = pont(ceil(jj));
            y10 = hpont/2;                                          % Coordinata y Nodo 10
            y11 = y10;                                              % Coordinata y Nodo 11

            mi_drawline(x11,y11,x10,y10);                           % Disegno base superiore ponticello
            mi_drawline(x11,-y11,x10,-y10);  
            
            x12 = x10 + 1.5 * racc_pont; y12 = y10+1.5*racc_pont;   % Punto 12 (serve solo per il disegno del raccordo: verrà cancellato)

            mi_drawline(x10,y10,x12,y12);                           % Disegno il raccordo di sx
            mi_drawline(x10,-y10,x12,-y12);

            mi_selectnode(x12,y12);                                 % Cancello i nodi di troppo (e così anche i tratti di segmento in eccesso)
            mi_selectnode(x12,-y12);
            mi_deleteselectednodes

            x13 = x11 - 1.5 * racc_pont; y13 = y10+1.5*racc_pont;   % Punto 13 (serve solo per il disegno del raccordo: verrà cancellato)

            mi_drawline(x11,y11,x13,y13);                           % Disegno il raccordo di dx
            mi_drawline(x11,-y11,x13,-y13);

            mi_selectnode(x13,y13);                                 % Cancello i nodi di troppo (e così anche i tratti di segmento in eccesso)
            mi_selectnode(x13,-y13);
            mi_deleteselectednodes
            
            mi_selectsegment(X32(jj),0);                             % Cancello l'archetto di confine sx
            mi_deleteselectedsegments;

            mi_selectsegment(X42(jj),0);                             % Cancello l'archetto di confine dx
            mi_deleteselectedsegments;

        else                                                         % Se invece i punti 10 e 11 si incrociano (ovvero raccordo troppo grande rispetto)
            % alla larghezza barriera), non faccio il ponticello e disegno solo una linea di
            mi_drawline(x11-racc_pont,0,x10+racc_pont,0);           % separazione. Questo vale anche nel caso in cui lo spessore impostato per il
            hpont = 0;                                              % ponticello sia troppo piccolo.

        end
        % DISEGNO DELLE BARRIERE DI FLUSSO -> Assegno magneti

        if (jj==1)
            x42=x2;
            x32=x1;
        x = x42+(x32-x42)/2;
        y = hpont;
        if hpont == 0
            y = 0.3*pont0;                                                                 % Eccezione: quando non c'è il ponticello radiale
        end
        else
           x = x42+(x32-x42)/2;
        y = hpont;
        if hpont == 0
            y = 0.3*pont0;                                                                 % Eccezione: quando non c'è il ponticello radiale
        end  
        end
        
        betac = calc_apertura_cerchio(alphacc(jj),rcc(jj),x0);
     
        mi_addblocklabel(x,y); mi_selectlabel(x,y);
      % mi_setblockprop(magnet, 1, 0, 'None', -beta(ceil(jj/2))/2 + 180, 0, 1);            % Nota: direzione di magnetizzazione
       mi_setblockprop(magnet, 1, 0, 'None', -135, group0, 1);
       %        mi_setblockprop(magnet, 1, 0, 'None',0, group0, 1);            % Nota: direzione di magnetizzazione
        mi_clearselected;

        mi_addblocklabel(x,-y); mi_selectlabel(x,-y);
        mi_setblockprop(magnet, 1, 0, 'None', -135, group0, 1);
       % mi_setblockprop(magnet, 1, 0, 'None',0, group0, 1);
        mi_clearselected;
        
        %% Set istruzioni per magneti radiali:
        % Si dividono i magneti in 3 segmenti in modo da attribuire la
        % magnetizzazione perpendicolare al segmento in questione:
        %% Decommentare le seguenti righe nel caso di magnetizzazione
        %% completa dei layer:
        if (jj~=1)
            mi_drawline(x32,y32,x42,y32);
            mi_drawline(x32,-y32,x42,-y32);
            mi_addblocklabel(x32,y42); mi_selectlabel(x32,y42);
            mi_setblockprop(magnet, 1, 0, 'None', 180, group0, 1);
            mi_clearselected;
            mi_addblocklabel(x32,-y42); mi_selectlabel(x32,-y42);
            mi_setblockprop(magnet, 1, 0, 'None', -90, group0, 1);
            mi_clearselected;

        end
        %% Nel caso in cui si usino i tegoli (commentare le seguenti righe
        %% per magnetizzazione completa layer e decommentare 382 to 392)
%         if (jj==1)
%             mi_drawline(X32(jj),Y32(jj),X42(jj),Y42(jj));
%             mi_drawline(X32(jj),-Y32(jj),X42(jj),-Y42(jj));
%             xa1=X32(jj)-(X32(jj)-X42(jj))/2;
%             ya1=Y32(jj)+hc(jj)/4;
%             mi_addblocklabel(xa1,ya1);mi_selectlabel(xa1,ya1);
%             mi_setblockprop('Air', 1, 0, 'None', 0, group0, 1);
%             mi_clearselected;
%             mi_addblocklabel(xa1,-ya1); mi_selectlabel(xa1,-ya1);
%             mi_setblockprop('Air', 1, 0, 'None', 0, group0, 1);
%             mi_clearselected;
%         elseif (jj==2)
% 
%             mi_drawline(x32,y32,x42,y32);
%             mi_drawline(x32,-y32,x42,-y32);
%             
%             nax_est=X3(jj);
%             nay_est=Y3(jj)+hc(jj)/2*cos(pi/4);
%             mi_addblocklabel(nax_est,nay_est);
%             mi_selectlabel(nax_est,nay_est);
%             mi_setblockprop('Air', 1, 0, 'None', 0, group0, 1);
%             mi_clearselected;
%             mi_addblocklabel(nax_est,-nay_est);
%             mi_selectlabel(nax_est,-nay_est);
%             mi_setblockprop('Air', 1, 0, 'None', 0, group0, 1);
%             mi_clearselected;
% 
%             
%         else
%             mi_drawline(x32,y32,x42,y32);
%             mi_drawline(x32,-y32,x42,-y32);
%             mi_addblocklabel(x32,y42); mi_selectlabel(x32,y42);
%             mi_setblockprop(magnet, 1, 0, 'None', 180, group0, 1);
%             mi_clearselected;
%             mi_addblocklabel(x32,-y42); mi_selectlabel(x32,-y42);
%             mi_setblockprop(magnet, 1, 0, 'None', -90, group0, 1);
%             mi_clearselected;
% 
%             %keyboard
%             mi_drawline(X3(jj),Y3(jj),X4(jj),Y4(jj))
%             mi_drawline(X3(jj),-Y3(jj),X4(jj),-Y4(jj))
%             
%             mi_drawline(X3(jj),Y3(jj),X4(jj),Y4(jj))
%             mi_drawline(X3(jj),-Y3(jj),X4(jj),-Y4(jj))
%             % coordinate tegoli radiali
%             nlayx=-hc(jj)*cos(pi/4)+X32(jj);
%             nlayy=hc(jj)*sin(pi/4)+Y32(jj);
%             mi_addnode(nlayx,nlayy);
%             mi_drawline(X32(jj),Y32(jj),nlayx,nlayy);
%             mi_drawline(X32(jj),-Y32(jj),nlayx,-nlayy);
%             %coordinate aria tra tegoli diritti e tegoli radiali
%             nax=-hc(jj)/2*cos(pi/4)+X32(jj);
%             nay=hc(jj)/4*sin(pi/4)+Y32(jj);
%             mi_addblocklabel(nax,nay);
%             mi_selectlabel(nax,nay);
%             mi_setblockprop('Air', 1, 0, 'None', 0, group0, 1);
%             mi_clearselected;
%             mi_addblocklabel(nax,-nay);
%             mi_selectlabel(nax,-nay);
%             mi_setblockprop('Air', 1, 0, 'None', 0, group0, 1);
%             mi_clearselected;
%             % aria vicino ai ponticelli nei layer radiali
%             nax_est=X3(jj);
%             nay_est=Y3(jj)+hc(jj)/2*cos(pi/4);
%             mi_addblocklabel(nax_est,nay_est);
%              mi_selectlabel(nax_est,nay_est);
%             mi_setblockprop('Air', 1, 0, 'None', 0, group0, 1);
%             mi_clearselected;
%             mi_addblocklabel(nax_est,-nay_est);
%             mi_selectlabel(nax_est,-nay_est);
%             mi_setblockprop('Air', 1, 0, 'None', 0, group0, 1);
%             mi_clearselected;
%         end
        %
%         tegoli_in_SmCo
        
        %%
   i=i+1;
   %Hfe(i)=hfe2;
   
end

%% (NOTA Matteo) RICORDA: non hai indicato in precedenza betac per la direzione dei
%% magneti, posto 0 in riga 356 e 360
%% magneti
% mi_addblocklabel(x,y); mi_selectlabel(x,y);
%       mi_setblockprop(magnet, 1, 0, 'None', -beta(ceil(jj/2))/2 + 180, 0, 1);            % Nota: direzione di magnetizzazione
%        mi_setblockprop(magnet, 1, 0, 'None', -betac/2 * 180/pi + 180, group0, 1);            % Nota: direzione di magnetizzazione
%        mi_clearselected;

%        mi_addblocklabel(x,-y); mi_selectlabel(x,-y);
%        mi_setblockprop(magnet, 1, 0, 'None', betac/2 * 180/pi + 180, group0, 1);
%        mi_clearselected;
%%

mi_clearselected;

%% Parametri che servono per valutare Pfe
% baricentro guide di flusso
r_fe = [r_all x0-Ar];
r_fe = r_fe(2:end);
% spessore ferri
hf = r_fe(2:2:end) - r_fe(1:2:end);
% coordinata angolare mediana ferri
alpha_f = [alpha 90/p];
alpha_f = mean([alpha_f(2:end);alpha_f(1:end-1)]);
% raggio mediano ferri (cerchio centrato in x0,0)

r_fe = mean([r_fe(1:2:end);r_fe(2:2:end)]);
% apertura angolare settore cerchio r_fe
beta_f = 180/pi * calc_apertura_cerchio(pi/180*alpha_f,xr,x0);
% punti centrali guide - nei pressi del punto 4
R4 = abs(X4 + 1i * Y4);
[x_fe,y_fe] = calc_intersezione_cerchi(R4,r_fe,x0);
%% calcolo dei punti centro guida sull'asse di simmetria (new geometry, Matteo 28/11/2011):
hf(:,nlay)=Xbar(nlay)-Ar;  %layer vicino albero modifica per beccare centro del segmento di ferro
r_fe=Xbar-hf/2;



%% COMPLETO IL DISEGNO (ARCO ESTERNO E TRAFERRO_STRATO1) E LE ASSEGNAZIONI (MATERIALI E CONDIZIONI AL CONTORNO)


% ARCO AL TRAFERRO

x = xr * cos(90/p * pi/180);
y = xr * sin(90/p * pi/180);

mi_drawarc(x,-y,x,y,180/p,1);

%keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % assegno magneti e ferro
% for jj = 1:nlay
%     if gammam(jj)>0
%         keyboard
%         x = sqrt(2) * xr - r(jj);
%         y = pont(jj);
%         % eccezione per magneti molto piccoli
%         if y > y7
%             y = mean([pont(jj)/2 y7]);
%         end
%         mi_addblocklabel(x,y); mi_selectlabel(x,y);
%         mi_setblockprop(magnet, 1, 0, 'None', -gammam(jj)/2 + 180, 0, 1);
%         mi_clearselected;
%         mi_addblocklabel(x,-y); mi_selectlabel(x,-y);
%         mi_setblockprop(magnet, 1, 0, 'None', gammam(jj)/2 + 180, 0, 1);
%         mi_clearselected;
%     end
% end
% mi_clearselected;
% assegno ferro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%keyboard
%% ASSEGNO FERRO
xf = xr - 0.5*pont0; yf = 0;
[xy_ferro(1),xy_ferro(2)] = rot_point(xf,yf,(90/p)*pi/180);
mi_addblocklabel(xf,yf);
mi_selectlabel(xf,yf);
mi_setblockprop(steel, 1, 0, 'None', 0, group0, 1);

mi_selectgroup(group0); mi_moverotate(0,0,90/p);                     % Ruoto nuovamente la figura

%keyboard
%% DISEGNO LO STRATO 1 DEL TRAFERRO (VA DA (th_m0) A (th_m0+180/p))

x1 = xr; y1 = 0;
x2 = xr+1/3*g; y2 = 0;
mi_drawline(x1,y1,x2,y2);
mi_selectsegment(mean([x1 x2]),mean([y1 y2]));
mi_setsegmentprop('APg1', 0, 1, 0, group0);

%keyboard

[x1,y1] = rot_point(x1,y1,180/p*pi/180);
[x2,y2] = rot_point(x2,y2,180/p*pi/180);
mi_drawline(x1,y1,x2,y2);
mi_selectsegment(mean([x1 x2]),mean([y1 y2]));
mi_setsegmentprop('APg1', 0, 1, 0, group0);

%keyboard

[x1,y1] = rot_point(xr+1/6*g,0,90/p*pi/180);
mi_addblocklabel(x1,y1);
mi_selectlabel(x1,y1);
mi_setblockprop('Air', 1, 0, 'None', 0, group0, 1);
mi_clearselected;
%keyboard
% NOTA BENE: Non disegno l'arco di circonferenza di raggio (xr+1/3*g) -> Vedi disegna_bordo_mobile

%% variabili da portare fuori
geo.pont = pont;
geo.r_all = r_all;

geo.hf = hf;
geo.beta_f = beta_f;
geo.r_fe = r_fe;
geo.dVol_Fe = dVol_Fe;
geo.x_fe = x_fe;
geo.y_fe = y_fe;
geo.Xbar = Xbar;
geo.X42 = X42;
geo.X32 = X32;
geo.Y42 = Y42;
geo.Y32 = Y32;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mi_saveas(filename) % saves the file with name ’filename’.
% mi_zoomnatural
% pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mi_selectgroup(group0); mi_setgroup(2);

keyboard

