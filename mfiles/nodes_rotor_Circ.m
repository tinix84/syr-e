% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

%% DISEGNO DEL ROTORE (GRUPPO 0)
%% modifica 15 novembre - gp
% - corretto angolo di magnetizzazione errato
% - aggiunto albero per evitare saturazione spider (irrealistica)
%% modifica 08 gennaio 2010 - gp
% - porto fuori geo.r_all (servono per Pfe)
% - calcolo e porto fuori geo.r_fe,geo.hf (servono per Pfe)

function [geo,temp] = nodes_rotor_Circ(geo) %,mat,fem)

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

dalpha = geo.dalpha;            % Angoli dalpha
% Eval alpha
alpha = integral_fx(1:length(dalpha),dalpha);
geo.alpha=alpha;
racc_pont = geo.racc_pont;      % racc_pont=1*pont0 <- per i ponticelli radiali.
ang_pont0 = geo.ang_pont0;      % Ampiezza dell'angolo (in gradi) da spazzare con  raggio xr in modo da ottenre un arco lungo pont0

nmax = geo.nmax;                % Velocità max (rpm) per la valutazione della sollecitazione centrifuga più gravosa (-> ponticelli)
error_mex=zeros(1,nlay);

%% determination of air thickness and check the feasibility of the geometry
geo = calcHcCheckGeoControl(geo);
hc=geo.hc;

%% CONTORNO DEL ROTORE: INDIVIDUAZIONE DELLE COORDINATE PRINCIPALI; DISEGNO
%% DELL'ALBERO E DELLE LINEE DI INZIO E FINE POLO.
% DEFINIZIONE DELLE COORDINATE CHE INDIVIDUANO IL CONTORNO PRINCIPALE DEL ROTORE VERSO IL TRAFERRO

x1 = xr;        y1 = 0;                             % Punto1
[x2,y2] = rot_point(x1,y1,180/p*pi/180);            % Punto2
x3 = xr-pont0;  y3 = 0;                             % Punto3
[x4,y4] = rot_point(x3,y3,180/p*pi/180);            % Punto4

x1 = Ar; y1 = 0;                                    % Punto1
[x2,y2] = rot_point(x1,y1,180/p*pi/180);            % Punto2

[x_shaft_mat,y_shaft_mat]=rot_point(x1/2,y1,90/p*pi/180);

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
XBanqdx=x0-r_all(1:2:end);
XBanqsx=x0-r_all(2:2:end);

%% CALCOLO PONTICELLI IN BASE ALLA CENTRIFUGA
% 1. volume guide di flusso
% 2. somma dei volumi e calcolo delle masse sospese su ciascun pont
% 3. calcolo dello spessore pont che corrisponde alla sollecitazione
% centrifuga voluta (sigma_max = 180 N/mm2)

% Volume singolo settore di ferro
% (i settori di ferro in cui si divide il rotore sono separati dalle circonferenze di raggio r centrate in (x0,0))
dVol_Fe = [];
rG = [];                        % Distanza da (0,0) dei baricentri dei volumi appesi ai ponticelli

for jj = 1:2:length(r_all)
    kk = ceil(jj/2);            % indice ponticello
    
    if kk == 1
        
        dVol_Fe(1) = l * (r(1)^2*(beta(1)*pi/180-0.5*sin(2*beta(1)*pi/180)) + xr^2*(alpha(1)*pi/180-0.5*sin(2*alpha(1)*pi/180)));
        
    else
        
        dVol_Fe(kk) = l * (r(kk)^2-r(kk-1)^2)*0.5*(beta(kk)*pi/180+beta(kk-1)*pi/180);
        
    end
    
end

Vol_Fe = integral_fx(1:length(dVol_Fe),dVol_Fe);      % volumi appesi ai ponticelli
M_Fe = Vol_Fe * 1e-9 * 7800 ;                        % massa ferro appeso ai ponticelli

rG =  0.5 * (xr + (x0 - r));
F_centrifuga = M_Fe .* rG/1000 *  (nmax * pi/30)^2;
sigma_max = 180;                                     % N/mm2 - snervamento lamierino

pont = F_centrifuga/(sigma_max * l);                 % mm

pont(pont < pont0) = 0;                              % NOTA BENE: Elimino i ponticelli troppo sottili

% DISEGNO DELLE BARRIERE DI FLUSSO (ASSEGNO MAGNETI E ARIA).
X4 = zeros(1,nlay); Y4 = X4;
X3 = zeros(1,nlay); Y3 = X4;

for jj = 1:2:length(r_all)                                                          % Le barriere da disegnare sono nlay. r_all ha invece dimensione
    % 2*nlay. Sarà quindi la quantità ceil(jj/2)a indicare, man mano,
    % il numero progressivo della barriera che si sta disegnando.
    % (Nota:il comando ceil arrotonda infatti il numero tra parentesi
    % all'intero maggiore o uguale più vicino).
    
    % DISEGNO DELLE BARRIERE DI FLUSSO -> Definizione del centro del semicerchio per il disegno della punta.
    
    if hc(ceil(jj/2))>0
        
        [xc,yc] = calc_intersezione_cerchi(xr-pont0-hc(ceil(jj/2))/2, r(ceil(jj/2)), x0);
        
        % DISEGNO DELLE BARRIERE DI FLUSSO -> Aggiunta dei nodi 3 (-3) e 4 (-4) (Bordi inferiore e superiore della barriera_inizio)
        
        % 2014/02/25 MG no intersaction between circle mean that the
        % barrier is unfeseable and is corresponding drawn like a circle
        if (not(isreal(xc))||not(isreal(yc)))
            xc=xr-pont0-hc(ceil(jj/2))/2;
            yc=0;
            x3=xc-hc(ceil(jj/2))/2; y3=yc;
            x4=xc+hc(ceil(jj/2))/2; y4=yc;
            error_mex(ceil(jj/2))=1;    % barrier not drawn;
        else
        [alphac,rc] = cart2pol(xc,yc);
        betac = calc_apertura_cerchio(alphac,rc,x0);
        [x3,y3] = calc_punto_magnete(r_all(jj), betac, x0);                     % Nodo 3 (Bordo inferiore della barriera_inizio)
        
        [x4,y4] = calc_punto_magnete(r_all(jj+1), betac, x0);                   % Nodo 4 (Bordo superiore della barriera_inizio)
        end    
        % 2013/07/06 MG si memorizzano in un vettore le coordinate dei
        % centri barriera
        XcBan(ceil(jj/2))=xc;   YcBan(ceil(jj/2))=yc;
        % memorizzo il punto 4 per il calcolo di Pfe
        X4(ceil(jj/2)) = x4;
        Y4(ceil(jj/2)) = y4;
        X3(ceil(jj/2))=x3;
        Y3(ceil(jj/2))=y3;
        
        if y3 >0
            
            [alpha3,r3] = cart2pol(x3,y3);
            beta3 = calc_apertura_cerchio(alpha3,r3,x0);                            % beta3=apertura angolare del punto (x3,y3) rispetto a (x0,0)
%             mi_drawarc(x3,y3,x3,-y3,2*beta3*180/pi,10);                             % Disegno del bordo inferiore
            
        end
        
        if y4 >0
            
            [alpha4,r4] = cart2pol(x4,y4);
            beta4 = calc_apertura_cerchio(alpha4,r4,x0);
%             mi_drawarc(x4,y4,x4,-y4,2*beta4*180/pi,10);                             % Disegno del bordo superiore
            
        end
        
        % DISEGNO DEI PONTICELLI RADIALI
        % keyboard
        x10 = x0 - r_all(jj) - racc_pont;                           % Coordinata x Nodo 10 ((Ponticello radiale_vertice in alto a dx)=(Punto
        % in basso a dx del raccordo))
        
        x11 = x0 - r_all(jj+1) + racc_pont;                         % Coordinata x Nodo 11 ((Ponticello radiale_vertice in alto a sx)=(Punto
        % in basso a sx del raccordo))
        
        if (x11 < x10) && (pont(ceil(jj/2))>0)                      % Se i due punti non si incrociano (raccordo non troppo grande rispetto alla
            % larghezza della barriera), e se lo spessore estratto dall'algoritmo per il
            % ponticello non è troppo piccolo, allora procedo al disegno del ponticello
            
            lpont = x10 - x11;                                      % lpont=lunghezza ponticello=larghezza barriera-2*racc_pont
            hpont = pont(ceil(jj/2));
            y10 = hpont/2;                                          % Coordinata y Nodo 10
            y11 = y10;                                              % Coordinata y Nodo 11
            
            XpontRadDx(ceil(jj/2))=x10; YpontRadDx(ceil(jj/2))=y10;
            XpontRadSx(ceil(jj/2))=x11; YpontRadSx(ceil(jj/2))=y10;
            
            x12 = x10 + 1.5 * racc_pont; y12 = y10+1.5*racc_pont;   % Punto 12 (serve solo per il disegno del raccordo: verrà cancellato)
            XpontRadBarDx(ceil(jj/2))=XBanqdx(ceil(jj/2)); YpontRadBarDx(ceil(jj/2))=y10+(XpontRadBarDx(ceil(jj/2))-x10);
                       
            mi_selectnode(x12,y12);                                 % Cancello i nodi di troppo (e così anche i tratti di segmento in eccesso)
            mi_selectnode(x12,-y12);
            mi_deleteselectednodes
            
            x13 = x11 - 1.5 * racc_pont; y13 = y10+1.5*racc_pont;   % Punto 13 (serve solo per il disegno del raccordo: verrà cancellato)
            XpontRadBarSx(ceil(jj/2))=XBanqsx(ceil(jj/2)); YpontRadBarSx(ceil(jj/2))=y10+(x11-XpontRadBarSx(ceil(jj/2)));
            
        else                                                         % Se invece i punti 10 e 11 si incrociano (ovvero raccordo troppo grande rispetto)
            % alla larghezza barriera), non faccio il ponticello e disegno solo una linea di
            hpont = 0;                                              
            XpontRadBarSx(ceil(jj/2))=XBanqsx(ceil(jj/2));
            YpontRadBarSx(ceil(jj/2))=0;
            XpontRadBarDx(ceil(jj/2))=XBanqdx(ceil(jj/2));
            YpontRadBarDx(ceil(jj/2))=0;
            XpontRadDx(ceil(jj/2))=NaN;
            YpontRadDx(ceil(jj/2))=0;
            XpontRadSx(ceil(jj/2))=NaN;
            YpontRadSx(ceil(jj/2))=0;
            
        end
        
        
        if ceil(jj/2) == 1;
            x = x0 - r(ceil(jj/2));
            y = hpont(1) * 0.6;
            if y == 0
                y = pont0;
            end
        else
            [xtemp,ytemp] = rot_point(-r(ceil(jj/2)),0,-0.5*beta(ceil(jj/2))*pi/180);
            x = x0 + xtemp;
            y = ytemp;
        end
    end
end
% mi_clearselected;
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

%% COMPLETO IL DISEGNO (ARCO ESTERNO E TRAFERRO_STRATO1) E LE ASSEGNAZIONI (MATERIALI E CONDIZIONI AL CONTORNO)
% ARCO AL TRAFERRO

x = xr * cos(90/p * pi/180);
y = xr * sin(90/p * pi/180);

% ASSEGNO FERRO
xf = xr - 0.5*pont0; yf = 0;
[xy_ferro(1),xy_ferro(2)] = rot_point(xf,yf,(90/p)*pi/180);

% DISEGNO LO STRATO 1 DEL TRAFERRO (VA DA (th_m0) A (th_m0+180/p))

x1 = xr; y1 = 0;
x2 = xr+1/3*g; y2 = 0;

[x1,y1] = rot_point(x1,y1,180/p*pi/180);
[x2,y2] = rot_point(x2,y2,180/p*pi/180);
% 
[x1,y1] = rot_point(xr+1/6*g,0,90/p*pi/180);

% %% variables to bring out
geo.pont = pont;
% geo.r_all = r_all;

geo.hf = hf;
% geo.beta_f = beta_f;
% geo.r_fe = r_fe;
% geo.dVol_Fe = dVol_Fe;
% geo.x_fe = x_fe;
% geo.y_fe = y_fe;

temp.XpontRadDx=XpontRadDx;
temp.YpontRadDx=YpontRadDx;
temp.XpontRadSx=XpontRadSx;
temp.YpontRadSx=YpontRadSx;
temp.XpontRadBarDx=XpontRadBarDx;
temp.XpontRadBarSx=XpontRadBarSx;
temp.YpontRadBarDx=YpontRadBarDx;
temp.YpontRadBarSx=YpontRadBarSx;

temp.XBanqdx=XBanqdx;
temp.XBanqsx=XBanqsx;
temp.XcBan=XcBan;
temp.YcBan=YcBan;
temp.X3=X3;
temp.Y3=Y3;

temp.X4=X4;
temp.Y4=Y4;
temp.error_mex=error_mex;

end
