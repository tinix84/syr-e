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

function [geo,mat,temp]=nodes_rotor_Circ_dx(geo,mat)

r = geo.r;                      % Raggio del rotore al traferro
x0 = geo.x0;                    % Centro fittizio
rshaft = geo.Ar;                % Raggio albero
Ar=geo.Ar;
l = geo.l;                      % Lunghezza pacco
g = geo.g;                      % Traferro
pont0 = geo.pont0;              % minimum mechanical tolerance
pontT = geo.pontT;              % Airgap ribs [mm]

p = geo.p;                      % Paia poli
nlay = geo.nlay;                % N° layers

dalpha = geo.dalpha;            % Angoli dalpha
% Eval alpha
alpha = cumsum(dalpha);
dx=geo.dx;

racc_pont = geo.racc_pont;      % racc_pont=1*pont0 <- per i ponticelli radiali.
ang_pont0 = geo.ang_pont0;      % Ampiezza dell'angolo (in gradi) da spazzare con  raggio r in modo da ottenre un arco lungo pont0

nmax = geo.nmax;                % Velocità max (rpm) per la valutazione della sollecitazione centrifuga più gravosa (-> ponticelli)
hfe_min=geo.hfe_min;
sigma_max = mat.Rotor.sigma_max;    % snervamento materiale [MPa]
rhoFE = mat.Rotor.kgm3;             % densità del ferro di rotore [kg/m3]
rhoPM = mat.LayerMag.kgm3;          % densità magneti [kg/m3]

XcRibTraf1=zeros(1,nlay);
XcRibTraf2=zeros(1,nlay);
YcRibTraf1=zeros(1,nlay);
YcRibTraf2=zeros(1,nlay);
xxD1k=zeros(1,nlay);
yyD1k=zeros(1,nlay);
xxD2k=zeros(1,nlay);
yyD2k=zeros(1,nlay);

% CENTRO FITTIZIO DI COORDINATE (x0,0)
beta = 180/pi * calc_apertura_cerchio(pi/180*alpha,r,x0);                  % La funzione calc_apertura_cerchio riceve in input le coordinate polari,
% (r, alpha) (alpha in rad), di un punto generico. Queste sono calcolate
% rispetto al centro (0,0). In output la funzione restituisce l'apertura
% angolare (in rad) dello stesso punto rispetto al centro preso come
% riferimento (ha coordinate:(x0,0)).
% I punti di cui, in questo caso, si calcolano le aperture angolari rispetto
% al centro di riferimento sono i punti mediani delle barriere, presi in
% corrispondenza del traferro.

rbeta = (x0 - r * cos(alpha*pi/180))./(cos(beta*pi/180));   % Di questi stessi punti, si calcolano anche le distanze dal centro (x0,0)
% e le si memorizzano nel vettore rbeta.

[xpont,ypont] = calc_intersezione_cerchi(r-pontT, rbeta, x0);

LowDimBarrier=zeros(1,nlay);
for ii=1:nlay
    if (not(isreal(xpont(ii)))||not(isreal(ypont(ii))))
        xpont(ii)=(r-2*pontT).*cos(alpha(ii)*pi/180);
        ypont(ii)=(r-2*pontT).*sin(alpha(ii)*pi/180);
        LowDimBarrier(ii)=1;
    end
end

rpont_x0=sqrt(ypont.^2+(x0-xpont).^2);
[alphapont,rpont] = cart2pol(xpont,ypont);
Bx0=x0-(rpont_x0);

% Determination of air thickness and check the feasibility of the geometry
geo.Bx0=Bx0; % Initialization of central non-moved line of the flux barrier
if geo.delta_FBS==0
    geo = calcHcCheckGeoControlwDx(geo);
    B1k=geo.B1k;
    B2k=geo.B2k;
else
    B1k=Bx0-geo.hc/2+geo.dx.*geo.hc/2;
    B2k=Bx0+geo.hc/2+geo.dx.*geo.hc/2;
    geo.B1k=B1k;
    geo.B2k=B2k;
end

hc=B2k-B1k;

ptmp=find(abs(hc)<pont0/2);
B2k(ptmp)=min([B1k(ptmp),B2k(ptmp)])+pont0/2;

% Intersezione circonferenze punti al traferro:
[xTraf2,yTraf2] = calc_intersezione_cerchi(r-pontT, x0-B2k, x0);
[xTraf1,yTraf1] = calc_intersezione_cerchi(r-pontT, x0-B1k, x0);

for ii=1:nlay
    if (not(isreal(xTraf1(ii)))||not(isreal(yTraf1(ii)))||not(isreal(xTraf2(ii)))||not(isreal(yTraf2(ii)))||hc(ii)<0)
        xpont(ii)=(r-2*pontT).*cos(alpha(ii)*pi/180);
        ypont(ii)=(r-2*pontT).*sin(alpha(ii)*pi/180);
        LowDimBarrier(ii)=1;
    end
end

% Barriers tips nodes (points 1 and 2 + centers)
for ii=1:nlay
    if LowDimBarrier(ii)==1
        xxD1k(ii)=B1k(ii);
        yyD1k(ii)=0;
        xxD2k(ii)=B2k(ii);
        yyD2k(ii)=0;
    else
        [xt,yt,xc,yc,rc]=cir_tg_2cir(xpont(ii),ypont(ii),r-pontT(ii),x0,0,x0-B1k(ii));
        xxD1k(ii)=xt;      % D1: end of outer semi arc
        yyD1k(ii)=yt;
        XcRibTraf1(ii)=xc; % center of semi arc 1
        YcRibTraf1(ii)=yc;
        [xt,yt,xc,yc,rc]=cir_tg_2cir(xpont(ii),ypont(ii),r-pontT(ii),x0,0,x0-B2k(ii));
        xxD2k(ii)=xt;      % D2: end of inner semi arc
        yyD2k(ii)=yt;
        XcRibTraf2(ii)=xc; % center of semi arc 2 
        YcRibTraf2(ii)=yc;
    end
end

temp.B1k=B1k;
temp.B2k=B2k;
temp.Bx0=Bx0;

[temp,geo] = calc_ribs_rad_fun(geo,mat,temp);

% Points for radial ribs
XpontRadDx=temp.XpontRadDx;
YpontRadDx=temp.YpontRadDx;
XpontRadSx=temp.XpontRadSx;
YpontRadSx=temp.YpontRadSx;
XpontRadBarDx=temp.XpontRadBarDx;
XpontRadBarSx=temp.XpontRadBarSx;
YpontRadBarDx=temp.YpontRadBarDx;
YpontRadBarSx=temp.YpontRadBarSx;

temp.xc=(xxD1k+xxD2k)/2;
temp.yc=(yyD1k+yyD2k)/2;

%% Determining  Magnet Area
YcBan = (yyD2k+yyD1k)/2;
XcBan = (xxD2k+xxD1k)/2;
eta1 = atand(yyD1k./(x0-xxD1k));
eta2 = atand(yyD2k./(x0-xxD2k));
eta = zeros(1,length(eta1));
for ii=1:length(eta1)
    if eta1(ii)<eta2(ii)
        eta(ii)=eta1(ii);
    else
        eta(ii)=eta2(ii);
    end
end

arc_area=pi*(eta/360).*[(x0-B1k).^2-(x0-B2k).^2]-([B2k-B1k].*YpontRadDx);   % surface of arc area for each barriers
pu_centerpost=([B2k-B1k].*YpontRadDx)./arc_area;                            %PU of center post area refer to arc area
Bar_fillfac=geo.BarFillFac;                                                 %barrier fill factor for real magnet

eta_min_dx=atand(YpontRadBarDx./(x0-XpontRadBarDx));
eta_min_sx=atand(YpontRadBarSx./(x0-XpontRadBarDx));
eta_min=zeros(1,length(eta_min_dx));
for ii=1:length(eta_min_dx)
    if eta_min_dx(ii)<eta_min_sx(ii)
        eta_min(ii)=eta_min_sx(ii);
    else
        eta_min(ii)=eta_min_dx(ii);
    end
end

%%%%%%%%%%%%%%%%%  detemining the ponits for real magnet
mag_angel=eta.*Bar_fillfac;
for ii=1:length(mag_angel)
    if mag_angel(ii)<eta_min(ii)
        mag_angel(ii)=eta_min(ii);
    end
end

R1 = x0-B2k;
R2 = x0-B1k;
%%%%%%%%%%%%determining coordinates for real magnet area
X5=x0-(cosd(mag_angel).*R1);
Y5=sind(mag_angel).*R1;
X6=x0-(cosd(mag_angel).*R2);
Y6=sind(mag_angel).*R2;

%%%%%%%%%%% Division of the PM in 2 part (to have a better approximation of
%%%%%%%%%%% the magnetization direction

half_angle = mag_angel/2;
for ii=1:length(half_angle)
    if half_angle(ii)<eta_min(ii)
        half_angle(ii)=eta_min(ii);
    end
end
xHalfDx = x0-(cosd(half_angle).*R1);
yHalfDx = sind(half_angle).*R1;
xHalfSx = x0-(cosd(half_angle).*R2);
yHalfSx = sind(half_angle).*R2;

%%%%%%%%%%%%%%
% determination of different magnet segment, central point and
% magnetization direction BarfillFac=0
xmag=[];
ymag=[];

if (sum(geo.BarFillFac)==0)
    xPMi = nan(size(Bx0));
    yPMi = nan(size(Bx0));
    xPMe = nan(size(Bx0));
    yPMe = nan(size(Bx0));
    
    for kk=1:nlay
        [a,b,c]=retta_per_2pti(XcBan(kk),YcBan(kk),x0,0);
        mOrto=-a/b/2;
        xmag=[xmag,cos(atan(mOrto))];
        ymag=[ymag,sin(atan(mOrto))];
    end
    mat.LayerMag.Br = [mat.LayerMag.Br mat.LayerMag.Br];    % replicates Br for correct block assignation
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % magnetization direction  BarfillFac~=0
    Xcc=(X5+X6)/2;
    Ycc=(Y5+Y6)/2;
    xPMi = x0-(x0-Bx0).*cosd(half_angle/2);
    yPMi = (x0-Bx0).*sind(half_angle/2);
    xPMe = x0-(x0-Bx0).*cosd(half_angle*3/2);
    yPMe = (x0-Bx0).*sind(half_angle*3/2);
    
    for kk=1:nlay
        xmag=[xmag,cosd(half_angle(kk)/2)];
        ymag=[ymag,-sind(half_angle(kk)/2)];
        xmag=[xmag,cosd(half_angle(kk)*3/2)];
        ymag=[ymag,-sind(half_angle(kk)*3/2)];
    end
    
    % ho 2 segmenti per ogni mezzo layer
    BrTemp=[];
    for ii=1:length(mat.LayerMag.Br)
        BrTemp=[BrTemp,mat.LayerMag.Br(ii) mat.LayerMag.Br(ii)];
    end
    mat.LayerMag.Br = BrTemp;
end

YcBan(YcBan==0)=hc(YcBan==0)/4;
temp.xmag=xmag;
temp.ymag=ymag;
temp.zmag=zeros(1,length(xmag));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results
temp.xpont=xpont;
temp.ypont=ypont;
temp.xxD1k=xxD1k;
temp.yyD1k=yyD1k;
temp.xxD2k=xxD2k;
temp.yyD2k=yyD2k;
temp.XcRibTraf1=XcRibTraf1;
temp.YcRibTraf1=YcRibTraf1;
temp.XcRibTraf2=XcRibTraf2;
temp.YcRibTraf2=YcRibTraf2;

% Points for BarFillFact
temp.X5=X5;
temp.Y5=Y5;
temp.X6=X6;
temp.Y6=Y6;

% Points for magnet segmentation
temp.xHalfDx = xHalfDx;
temp.yHalfDx = yHalfDx;
temp.xHalfSx = xHalfSx;
temp.yHalfSx = yHalfSx;

% Center points of the PM segments
temp.xPMi = xPMi;
temp.yPMi = yPMi;
temp.xPMe = xPMe;
temp.yPMe = yPMe;

temp.XcBan=XcBan;
temp.YcBan=YcBan;

% add control if the construction of barrier is too small, so barrier is
% constructed like diamond
temp.LowDimBarrier=LowDimBarrier;
hf=[r,B1k]-[B2k,Ar]; %calcolo dei Delta fi ferro di rotore
geo.hf = hf;
% geo.pont = pont;
geo.hc=hc;
geo.xpont=xpont;
geo.ypont=ypont;

% barrier transverse dimension (for permeance evaluation)
temp_r_arcs = (geo.r_all(1:2:end-1)+geo.r_all(2:2:end))/2;
temp_x_arc_ends = x0 - (xpont);
% equivalent sk = arc barrier, disregarding the end barrier radius
sk = temp_r_arcs.*acos(temp_x_arc_ends./temp_r_arcs);

geo.sk = sk;
geo.pbk = geo.sk ./ geo.hc;
geo.la = sum(geo.hc)/geo.r;
geo.lfe = sum(geo.hf)/geo.r;
geo.ly = (geo.R - (geo.r + geo.g + geo.lt))/geo.r;
geo.B1k=B1k;
geo.B2k=B2k;
geo.hc=hc;

