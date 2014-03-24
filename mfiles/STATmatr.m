%%
%% Stator
%%
% described by matrixes contatining 
% the master nodes
% the description of straight lines
% the description of arches
% the nodes where block labels must be placed

ns = geo.ns;
b = geo.b;
kt = geo.kt;
acs = geo.acs;
avv = geo.avv;
Nbob = geo.Nbob;
lt = geo.lt;
p = geo.p;
r = geo.r;
xr = geo.xr;
g = geo.g;
% Kbuccia = geo.kbuccia;
ttd = geo.ttd;
tta = geo.tta;
RaccordoFC = geo.SFR;

% mesh resolution
mesh=fem.res;
res_traf=fem.res_traf;

s.AS1=360/(ns*p)/2;             % angolo di mezzo passo cava
s.RSI=xr+g;                     % r traferro statore
LDente=kt*b*s.AS1*pi/180*s.RSI; % semilarghezza dente come metà arco
s.RSE=r;                        % r esterno statore
s.RS5=xr+g+lt;                  % r esterno cave
PassoCava=(s.RSI)*2*pi/(ns*p);  % passo cava di statore

s.C0=acs*PassoCava/2;
% s.C3=Kbuccia*b*(1-b*kt)*s.AS1*pi/180*s.RSI;
% s.DS2=b/2*PassoCava*acs;      % automatically evaluated...
s.DS2=ttd;                   % operator setup...


s.RS3=s.RSI+s.DS2;
[xtemp,ytemp] = rot_point(s.RS3,-LDente,s.AS1*pi/180);
s.DS1=sqrt((s.RS3-xtemp)^2+ytemp^2)*tan(tta*pi/180);
% s.DS1=s.C3-s.DS2;
s.C3=s.DS1+s.DS2;

s.RS1 = s.RSI+s.C3;            %tut. pag.3.31 (d0/c0=b/2)
[x9,y9] = rot_point(s.RS1,-LDente,s.AS1*pi/180);
s.RC1 = y9;

RC2o=s.RC1+(s.RS5-s.RS1)*tan(pi/180*s.AS1);     %semi_larghezza fondo cava SENZA arrotondamento
s.DS3=s.RS5*tan(s.AS1/2*pi/180);                %fondo cava

[x12,y12] = rot_point(s.RS5,-LDente,s.AS1*pi/180);
s.RS2=x12;     %fondo cava
s.RC2=y12;     %fondo cava

s.XS1=s.DS2+s.DS1+((s.RC1-s.DS1)^2-(s.DS1)^2)/(2*s.DS1);    %centro raccordo espansione dente_X
s.YS1=s.DS1;                                                %centro raccordo espansione dente_Y
s.XS2=s.RS2-s.RSI;
s.YS2=0;

s.RS4=(s.RS2+s.RS1)/2;                                      %separatore 1°- 2°strato
s.RC3=(s.RC1+s.RC2)/2;

[xRac,yRac,xtan1,ytan1,xtan2,ytan2]=cir_tg_2rette(s.RS5,0,s.RS5,s.DS3,s.RS1,s.RC1,s.RS4,s.RC3,RaccordoFC);
xRacFC=xRac(4); yRacFC=yRac(4);
xtgFC1=xtan1(4); ytgFC1=ytan1(4);
xtgFC2=xtan2(4); ytgFC2=ytan2(4);

s.filly1=[0; 0; s.DS3; s.RC2; s.RC3; s.RC1; s.DS1; s.DS1];
s.fillx=[s.RSI; s.RS5; s.RS5; s.RS2; s.RS4; s.RS1; s.RS3; s.RSI];
s.filly2=-s.filly1;

% Nodes
[x3,y3]=pol2cart(s.AS1*pi/180*acs,s.RSI);
[x4,y4]=pol2cart(s.AS1*pi/180,s.RSI);
[x7,y7]=pol2cart(s.AS1*pi/180,s.RS3);
[x10,y10]=pol2cart(s.AS1*pi/180,s.RS1);
[x13,y13]=pol2cart(s.AS1*pi/180,s.RS2);
[x16,y16]=pol2cart(s.AS1*pi/180,s.RSE);
[x19,y19]=pol2cart(s.AS1*pi/180,s.RS4);
[x20,y20]=pol2cart(s.AS1*pi/180,s.RS5);
xs1=s.RSI; ys1=0;
xs2=x3; ys2=y3;
xs5=x4; ys5=y4;
xs6=s.RSI+s.DS2; ys6=s.C0;
xs7=s.RS1; ys7=s.RC1;
xs8=xtgFC1; ys8=ytgFC1;
xs9=xtgFC2; ys9=ytgFC2;
xs10=s.RS4; ys10=s.RC3;
xs11=s.RS5; ys11=0;
xs12=s.RS1; ys12=0;
xs13=x10; ys13=y10;
xs14=s.RS4; ys14=0;
xs15=x4; ys15=y4;
xs16=x16; ys16=y16;
xs17=s.RSE; ys17=0;
xs89=xRacFC; ys89=yRacFC;

% Matrix describing lines and arches of half slot
CavaMat=[0 0 xs1 ys1 xs2 ys2 1;
    0 0 xs2 ys2 xs5 ys5 1;
    xs2 ys2 xs6 ys6 NaN NaN 0;
    xs6 ys6 xs7 ys7 NaN  NaN 0;
    %xs7 ys7 xs13 ys13 NaN  NaN 0;
    %xs13 ys13 xs15 ys15 NaN  NaN 0;
    xs15 ys15 xs16 ys16 NaN  NaN 0;
    %xs13 ys13 xs16 ys16 NaN  NaN 0;
    0 0 xs17 ys17 xs16 ys16  1;
    xs8 ys8 xs11 ys11 NaN  NaN 0;
    xs89 ys89 xs8 ys8 xs9 ys9 1;
    %         xs7 ys7 xs10 ys10 NaN  NaN 0 mesh;
    %         xs10 ys10 xs9 ys9 NaN  NaN 0 mesh;
    %         xs7 ys7 xs12 ys12 NaN  NaN 0 mesh;
    %         xs10 ys10 xs14 ys14 NaN  NaN 0 mesh
    ];

% lines dividing the slot copper in 2 or more layers (bundles)
Num_stat_strati=size(avv,1);
CavaStr=[];
for jk=1:Num_stat_strati
    xstr(jk)=s.RS1+(jk-1)*(s.RS5-s.RS1)/Num_stat_strati;
    ystr(jk)=s.RC1+(s.RC2-s.RC1)/(s.RS2-s.RS1)*(xstr(jk)-s.RS1);
    CavaStr=[CavaStr;[xstr(jk),0,xstr(jk),ystr(jk),NaN,NaN,0]];
    
end
CavaMat=[CavaMat;CavaStr];
CavaMat=[CavaMat;[xs7 ys7 xstr(1) ystr(1) NaN  NaN 0]];

for jk=1:Num_stat_strati-1
    CavaMat=[CavaMat;[xstr(jk),ystr(jk),xstr(jk+1),ystr(jk+1),NaN,NaN,0]];
    
end
CavaMat=[CavaMat;[xstr(end),ystr(end),xs9,ys9,NaN,NaN,0]];

% MIRROR of the half slot -> entire slot
CavaMatNeg=CavaMat;
CavaMatNeg(:,2)=-CavaMat(:,2);
CavaMatNeg(:,4)=-CavaMatNeg(:,4);
CavaMatNeg(:,6)=-CavaMat(:,6);
CavaMatNeg(:,7)=-CavaMat(:,7);

CavaMat=[CavaMatNeg;CavaMat];

% replicate/rotate one slot to obtain all slots
[nrig,ncol] = size(CavaMat);
MatRot=[];
for num_cave=1:(Qs-1)
    jk=1;
    CavaRot=[];
    while jk<=(ncol-2)
        [xrot,yrot]=rot_point(CavaMat(:,jk),CavaMat(:,jk+1),2*num_cave*s.AS1*pi/180);
        CavaRot=[CavaRot,xrot,yrot];
        jk=jk+2;
        
    end
    CavaRot=[CavaRot,CavaMat(:,ncol)];
    MatRot=[MatRot;CavaRot];
end
CavaMat=[CavaMat;MatRot];

% Materials codes
codAir = 2;
codFeStat = 4;
codCuStat = 3;

%% Block labels

% slot air
x_air_slot=(s.RSI+s.RS1)/2;
y_air_slot=0;
Air_slot=[x_air_slot,y_air_slot];
Temp=[];
nome_air_slot=[];
for num_cave=1:(Qs-1)
    
    [xtemp,ytemp]=rot_point(Air_slot(:,1),Air_slot(:,2),2*num_cave*s.AS1*pi/180);
    Temp=[Temp;xtemp,ytemp];
    
    %    nome_air_slot=[nome_air_slot;['slot_air_',num2str(num_cave)]];
    %    nome_air_slot_cell{num_cave}=nome_air_slot;
end

for num_cave=1:Qs
    
    %        nome_air_slot=[nome_air_slot;['slot_air_',num2str(num_cave)]];
    nome_air_slot{num_cave}={['slot_air_',num2str(num_cave)]};
    
end

Air_slot=[Air_slot;Temp];
[nrig,ncol]=size(Air_slot);
Air_slot=[Air_slot,codAir*ones(nrig,1),mesh*ones(nrig,1),0*ones(nrig,1)];
clear Temp xtemp ytemp

% tooth shoe
[xFe_esp_pos,yFe_esp_pos]=rot_point((s.RSI+s.RS1)/2,-LDente/2,s.AS1*pi/180);
xFe_esp_neg=xFe_esp_pos;
yFe_esp_neg=-yFe_esp_pos;
Fe_esp=[xFe_esp_pos,yFe_esp_pos;xFe_esp_neg,yFe_esp_neg];
Temp=[];
nome_Fe_esp=[];
for num_cave=1:(Qs-1)
    
    [xtemp,ytemp]=rot_point(Fe_esp(:,1),Fe_esp(:,2),2*num_cave*s.AS1*pi/180);
    Temp=[Temp;xtemp,ytemp];
    
end

Fe_esp=[Fe_esp;Temp];
[nrig,ncol]=size(Fe_esp);
Fe_esp=[Fe_esp,codFeStat*ones(nrig,1),mesh*ones(nrig,1),0*ones(nrig,1)];
clear Temp xtemp ytemp

for kk=1:size(Fe_esp,1)
    nome_Fe_esp{kk}={['Fe_esp_',num2str(kk)]};
end

% copper
% Cu_slot=[numero di cave x 2];
x_Cu_piu=[xstr(2:end),s.RS5];
x_Cu_min=xstr;
x_Cu_slot=(x_Cu_piu+x_Cu_min)/2;
y_Cu_slot=zeros(1,length(ystr));
Cu_slot=[x_Cu_slot',y_Cu_slot'];
Temp=[];
nome_Cu_slot=[];

for num_cave=1:(Qs-1)
    
    [xtemp,ytemp]=rot_point(x_Cu_slot',y_Cu_slot',2*num_cave*s.AS1*pi/180);
    Temp=[Temp;xtemp,ytemp];
    
end

ii=1;
for num_cave=1:Qs
    
    for num_str=1:size(avv,1)
        nome_Cu_slot{ii}={['slot_',num2str(num_str),'_',num2str(num_cave)]};
        ii=ii+1;
    end
    
end

Cu_slot=[Cu_slot;Temp];
[nrig,ncol]=size(Cu_slot);
Cu_slot=[Cu_slot,codCuStat*ones(nrig,1),mesh*ones(nrig,1),0*ones(nrig,1)];
clear Temp xtemp ytemp

% back iron
xFeYoke=(s.RS5+s.RSE)/2;
yFeYoke=0;

FeYoke=[xFeYoke,yFeYoke];
Temp=[];
nome_FeYoke=[];

for num_cave=1:(Qs-1)
    
    [xtemp,ytemp]=rot_point(xFeYoke,yFeYoke,2*num_cave*s.AS1*pi/180);
    Temp=[Temp;xtemp,ytemp];
    %     nome_FeYoke=[nome_FeYoke;['statore_',num2str(num_cave)]];
    
end

for num_cave=1:Qs
    nome_FeYoke{num_cave}={['statore_',num2str(num_cave)]};
end

FeYoke=[FeYoke;Temp];
[nrig,ncol]=size(FeYoke);
FeYoke=[FeYoke,codFeStat*ones(nrig,1),mesh*ones(nrig,1),0*ones(nrig,1)];
clear Temp xtemp ytemp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Boundary conditions

% FEMM-style codes
codBound_FluxTan  = 0;      % 0 flux tangential
codBound_periodic = 10;     % 10 odd or even periodicity

[xBoundRSE1,yBoundRSE1]=rot_point(s.RSE,0,s.AS1/2*pi/180);  % rotazione di 1/4 di passo cava
[xBoundRSE2,yBoundRSE2]=rot_point(s.RSE,0,-s.AS1/2*pi/180); % rotazione di 1/4 di passo cava
BoundRSE=[[xBoundRSE1,yBoundRSE1];[xBoundRSE2,yBoundRSE2]];
Temp=[];
for num_cave=1:Qs-1
    [xtemp,ytemp]=rot_point([xBoundRSE1;xBoundRSE2],[yBoundRSE1;yBoundRSE2],num_cave*2*s.AS1*pi/180);
    Temp=[Temp;xtemp,ytemp];
end

BoundRSE=[BoundRSE;Temp];
clear Temp xtemp ytemp;
BoundRSE=[BoundRSE,codBound_FluxTan*ones(size(BoundRSE,1),1)];

[xBoundLAT1,yBoundLAT1]=rot_point(mean([s.RSI,s.RSE]),0,-s.AS1*pi/180);
[xBoundLAT2,yBoundLAT2]=rot_point(mean([s.RSI,s.RSE]),0,(2*Qs-1)*s.AS1*pi/180);

BoundLAT=[xBoundLAT1,yBoundLAT1,codBound_periodic;xBoundLAT2,yBoundLAT2,codBound_periodic];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output matrixes:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% matrix 1: lines and arches
statore=CavaMat;

% matrix 2: block labels
BLKLABELSstat.xy=[Air_slot;Cu_slot;FeYoke];
BLKLABELSstat.names.air_slot=nome_air_slot';
BLKLABELSstat.names.Cu_slot=nome_Cu_slot';
BLKLABELSstat.names.FeYoke=nome_FeYoke';
BLKLABELSstat.names.legend={'air_slot','Cu_slot','FeYoke'};

% matrix 3: buondary conditions
BLKLABELSstat.boundary=[BoundLAT;BoundRSE];

% other data
STATOREdati.DatiVari.RSE=s.RSE;
STATOREdati.DatiVari.RSI=s.RSI;
STATOREdati.DatiVari.AngPC=s.AS1*2;
STATOREdati.DatiVari.Qs=Qs;


