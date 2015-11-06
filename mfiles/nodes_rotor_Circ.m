% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function [geo,temp] = nodes_rotor_Circ(geo) %,mat,fem)

group0 = 0;                      % di default il rotore è gruppo 0, alla fine dello script gli assegno il gruppo 2 (mod. 28 04 10)

% NOMENCLATURA: PROPRIETA' GEOMETRICHE E MATERIALI

r = geo.r;                      % Raggio del rotore al traferro
x0 = geo.x0;                    % Centro fittizio
Ar = geo.Ar;                    % Raggio albero
l = geo.l;                      % Lunghezza pacco
g = geo.g;                      % Traferro
pont0 = geo.pont0;              % Ponticelli al traferro (i ponticelli al traferro hanno lo spessore di un arco lungo pont0)

p = geo.p;                      % Paia poli
nlay = geo.nlay;                % N° layers

dalpha = geo.dalpha;            % Angoli dalpha
% Eval alpha
alpha = cumsum(dalpha);
racc_pont = geo.racc_pont;      % racc_pont=1*pont0 <- per i ponticelli radiali.
ang_pont0 = geo.ang_pont0;      % Ampiezza dell'angolo (in gradi) da spazzare con  raggio r in modo da ottenre un arco lungo pont0

nmax = geo.nmax;                % Velocità max (rpm) per la valutazione della sollecitazione centrifuga più gravosa (-> ponticelli)
error_mex=zeros(1,nlay);
sigma_max = geo.sigma_max;

% determination of air thickness and check the feasibility of the geometry
geo = calcHcCheckGeoControl(geo);
hc=geo.hc;

%% CONTORNO DEL ROTORE: INDIVIDUAZIONE DELLE COORDINATE PRINCIPALI; DISEGNO
%% DELL'ALBERO E DELLE LINEE DI INZIO E FINE POLO.
% DEFINIZIONE DELLE COORDINATE CHE INDIVIDUANO IL CONTORNO PRINCIPALE DEL ROTORE VERSO IL TRAFERRO

x1 = r;        y1 = 0;                             % Punto1
[x2,y2] = rot_point(x1,y1,180/p*pi/180);           % Punto2
x3 = r-pont0;  y3 = 0;                             % Punto3
[x4,y4] = rot_point(x3,y3,180/p*pi/180);           % Punto4

x1 = Ar; y1 = 0;                                   % Punto1
[x2,y2] = rot_point(x1,y1,180/p*pi/180);           % Punto2

[x_shaft_mat,y_shaft_mat]=rot_point(x1/2,y1,90/p*pi/180);

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
% rbeta = (x0 - r * cos(alpha*pi/180))./(cos(beta*pi/180));                      % Di questi stessi punti, si calcolano anche le distanze dal centro (x0,0)
B1k_temp = geo.B1k;
B2k_temp = geo.B2k;

r_all = [];                                                                 % Nel vettore r_all si memorizzano invece le distanze, sempre rispetto al
% Starting and ending point of the flux barrier...
for jj = 1:nlay                                                             % che individuano l'inizio e la fine delle nlay barriere.
    r_all = [r_all x0-B2k_temp(jj) x0-B1k_temp(jj)];
end
% median point of flux barrier...
rbeta=(r_all(1:2:end)+r_all(2:2:end))./2;
%% Posizione banane su asse q (convenzione assi VAGATI)
XBanqdx=x0-r_all(1:2:end);
XBanqsx=x0-r_all(2:2:end);

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
        
        [xc,yc] = calc_intersezione_cerchi(r-pont0-hc(ceil(jj/2))/2, rbeta(ceil(jj/2)), x0);
        %         [xp,yp]=calc_intersezione_cerchi(r-pont0, rbeta(ceil(jj/2)), x0);
        % DISEGNO DELLE BARRIERE DI FLUSSO -> Aggiunta dei nodi 3 (-3) e 4 (-4) (Bordi inferiore e superiore della barriera_inizio)
        
        % 2014/02/25 MG no intersaction between circle mean that the
        % barrier is unfeseable and is corresponding drawn like a circle
        if (not(isreal(xc))||not(isreal(yc))||yc<=eps)
            xc=r-pont0-hc(ceil(jj/2))/2;
            yc=0;
            x3=xc-hc(ceil(jj/2))/2; y3=yc;
            x4=xc+hc(ceil(jj/2))/2; y4=yc;
            error_mex(ceil(jj/2))=1;    % barrier not drawn;
        else
            [alphac,rc] = cart2pol(xc,yc);
            betac = calc_apertura_cerchio(alphac,rc,x0);
            %             [x3,y3,xc3,yc3,~]=cir_tg_2cir(xp,yp,r-pont0,x0,100*eps,r_all(jj)) ;
            [x3,y3] = calc_punto_magnete(r_all(jj), betac, x0);                     % Nodo 3 (Bordo inferiore della barriera_inizio)
            %             [x4,y4,xc4,yc4,~]=cir_tg_2cir(xp,yp,r-pont0,x0,100*eps,r_all(jj+1)) ;
            [x4,y4] = calc_punto_magnete(r_all(jj+1), betac, x0);                   % Nodo 4 (Bordo superiore della barriera_inizio)
            %             xc2=mean([x3,x4]); yc2=mean([y3,y4]);
        end
        % Barriers centers
        XcBan(ceil(jj/2))=xc;
        YcBan(ceil(jj/2))=yc;
        
        %         Xc3(ceil(jj/2))=xc3;   Yc3(ceil(jj/2))=yc3;
        %         Xc4(ceil(jj/2))=xc4;   Yc4(ceil(jj/2))=yc4;
        %         Xp(ceil(jj/2))=xp;   Yp(ceil(jj/2))=yp;
        % memorizzo il punto 4 per il calcolo di Pfe
        X4(ceil(jj/2)) = x4;
        Y4(ceil(jj/2)) = y4;
        X3(ceil(jj/2)) = x3;
        Y3(ceil(jj/2)) = y3;
        
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
    end
end

rTemp = rbeta;

calc_ribs_rad;

% determination of different magnet segment, central point and
% magnetization direction
xmag=[];
ymag=[];

for kk=1:nlay
    [a,b,c]=retta_per_2pti(XcBan(kk),YcBan(kk),x0,0);
    mOrto=-a/b/2;
    xmag=[xmag,cos(atan(mOrto))];
    ymag=[ymag,sin(atan(mOrto))];
        
end

geo.Br = [geo.Br geo.Br];    % replicates Br for correct block assignation

YcBan(YcBan==0)=hc(YcBan==0)/4;
temp.xmag=xmag;
temp.ymag=ymag;
temp.zmag=zeros(1,nlay);

% Parametri che servono per valutare Pfe
% baricentro guide di flusso
r_fe = [r_all x0-Ar];
r_fe = r_fe(2:end);
% spessore ferri
hf = r_fe(2:2:end) - r_fe(1:2:end);

% Output variables
geo.pont = pont;
geo.hf = zeros(1,nlay);
geo.hf = hf;

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

% barrier transverse dimension (for permeance evaluation)
temp_r_arcs = (r_all(1:2:end-1)+r_all(2:2:end))/2;
temp_x_arc_ends = x0 - (X3+X4)/2;
% equivalent sk = arc + barrier end radius weighted by 1/0.7822
sk = temp_r_arcs.*acos(temp_x_arc_ends./temp_r_arcs) + 1/0.7822 * hc/2;

geo.sk = sk;
geo.pbk = geo.sk ./ geo.hc;
geo.la = sum(geo.hc)/geo.r;
geo.lfe = sum(geo.hf)/geo.r;
geo.ly = (geo.R - (geo.r + geo.g + geo.lt))/geo.r;

temp.xc = temp.XcBan;
temp.yc = temp.YcBan;



