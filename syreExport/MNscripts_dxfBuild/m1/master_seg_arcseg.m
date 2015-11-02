%%%%%%%%%%%%%%%%%%%
% Parametri per il calcolo ed il disegno 
%               dello Statore
%%%%%%%%%%%%%%%%%%%%%
%
% cava trapezia definita secondo i parametri dei campisti
%

s.AS1=360/(ns*p)/2;             % angolo di mezzo passo cava
s.RSI=xr+g;                     % r traferro statore
LDente=kt*b*s.AS1*pi/180*s.RSI; % semilarghezza dente come metà arco
s.RSE=r;                        % r esterno statore
s.RS5=xr+g+lt;                  % r esterno cave
PassoCava=(s.RSI)*2*pi/(ns*p);  % passo cava di statore

s.DS1=acs*PassoCava/2;          % semi-apertura cava           
s.DS2=s.DS1;                    % provvisorio

% s.DS1=acs*PassoCava/2;          % semi-apertura cava           
% s.DS2=s.DS1/1.6;                    % provvisorio

s.RS3=s.RSI+s.DS2;

% old
% s.RS1 = s.RS3+s.DS2;          %tut. pag.3.31 (d0/c0=b/2)
% s.RC1 = (PassoCava-2*LDente)/2*s.RS1/s.RSI;

% corretto 24 08 09
s.RS1 = s.RS3+s.DS2;            %tut. pag.3.31 (d0/c0=b/2)
[x9,y9] = rot_point(s.RS1,-LDente,s.AS1*pi/180); 
s.RC1 = y9;

RC2o=s.RC1+(s.RS5-s.RS1)*tan(pi/180*s.AS1); %semi_larghezza fondo cava SENZA arrotondamento
s.DS3=s.RS5*tan(s.AS1/2*pi/180);                %fondo cava

% old
% s.RS2=s.RS5*cos(s.AS1*pi/180)+LDente*sin(s.AS1*pi/180);     %fondo cava
% s.RC2=s.RS5*sin(s.AS1*pi/180)-LDente*cos(s.AS1*pi/180);     %fondo cava

% corretto 24 08 09
[x12,y12] = rot_point(s.RS5,-LDente,s.AS1*pi/180);
s.RS2=x12;     %fondo cava
s.RC2=y12;     %fondo cava

s.XS1=s.DS2+s.DS1+((s.RC1-s.DS1)^2-(s.DS1)^2)/(2*s.DS1);    %centro raccordo espansione dente_X
s.YS1=s.DS1;                                                %centro raccordo espansione dente_Y      
s.XS2=s.RS2-s.RSI;
% s.YS2=s.DS3-((s.RC2-s.DS3)/2+raccordoD^2/(2*(s.RC2-s.DS3)));
s.YS2=0;

s.RS4=(s.RS2+s.RS1)/2;                                      %separatore 1°- 2°strato
s.RC3=(s.RC1+s.RC2)/2;

s.filly1=[0; 0; s.DS3; s.RC2; s.RC3; s.RC1; s.DS1; s.DS1];
s.fillx=[s.RSI; s.RS5; s.RS5; s.RS2; s.RS4; s.RS1; s.RS3; s.RSI];
s.filly2=-s.filly1;

% definizione della matrice dei nodi master
% statore;    % carica parametri statore
[x3,y3]=pol2cart(s.AS1*pi/180*acs,s.RSI);
[x4,y4]=pol2cart(s.AS1*pi/180,s.RSI);
[x7,y7]=pol2cart(s.AS1*pi/180,s.RS3);
[x10,y10]=pol2cart(s.AS1*pi/180,s.RS1);
[x13,y13]=pol2cart(s.AS1*pi/180,s.RS2);
[x16,y16]=pol2cart(s.AS1*pi/180,s.RSE);
[x19,y19]=pol2cart(s.AS1*pi/180,s.RS4);
[x20,y20]=pol2cart(s.AS1*pi/180,s.RS5);

mast=[];
mast=[0 0       %1
   s.RSI 0
   x3 y3
   x4 y4
   s.RSI+s.DS2 0    %5
   s.RSI+s.DS2 s.DS1
   x7 y7    %7
   s.RS1 0
   s.RS1 s.RC1
   x10 y10  %10
   s.RS2 0
   s.RS2 s.RC2
   x13 y13
   s.RS5 0
   s.RSE 0
   x16 y16
   s.RS4 0
   s.RS4 s.RC3  %18
   x19 y19
   x20 y20
   0 0
   0 0
   s.RS5 s.DS3];    %23
nummast=length(mast);

% segmenti

segm =  [3 6             % fianco dente interno
        6 9             % raccordo squadrato fondo cava
        9 18             % fondo cava verticale
        18 12
        12 23
        23 14
        8 9
        9 10
        17 18
        4 10
        10 16];

numseg = length(segm);

% archi
arcseg=[2 3 s.AS1*acs 0.5 0         % semiapertura cava (grid 0.5deg)
        3 4 s.AS1*(1-acs) 0.5 0     % mezzo dente (grid 0.5deg)
        15 16 s.AS1 s.AS1/5 0];

numarcseg=size(arcseg,1);