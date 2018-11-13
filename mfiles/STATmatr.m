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

function [geo,statore,BLKLABELSstat] = STATmatr(geo,fem)
%%
%% Stator
%%
% described by matrixes contatining 
% the master nodes
% the description of straight lines
% the description of arches
% the nodes where block labels must be placed

Qs = geo.Qs;
ns = geo.ns;
acs = geo.acs;
avv = geo.avv;
Nbob = geo.Nbob;
lt = geo.lt;
wt = geo.wt;
p = geo.p;
q = geo.q;
R = geo.R;
r = geo.r;
g = geo.g;
ttd = geo.ttd;
tta = geo.tta;
RaccordoFC = geo.SFR;

% mesh resolution
mesh=fem.res;
res_traf=fem.res_traf;

alpha_slot=2*pi/(ns*p);          % angolo di mezzo passo cava
RSI=r+g;                     % r traferro statore
RSE=R;                        % r esterno statore
RS5=r+g+lt;                  % r esterno cave
PassoCava=(RSI)*2*pi/(ns*p);  % passo cava di statore

% C0=acs*PassoCava/2;
% DS2=ttd;                   % operator setup...

% r eq for middle of slot computation
mr=tan(alpha_slot/2);
x0=RSI*cos(alpha_slot/2);
y0=RSI*sin(alpha_slot/2);
% tooth side computation
% design like a line parallel to r, case of trapezoidal slot r2: y=m2x+q2
% explicit form
m2=mr; q2= -wt/2*sqrt(1+mr^2);
% xx=linspace(0,100);
% y1=mr*xx;
% y2=m2*xx+q2;

[x1t,y1t]=intersezione_retta_circonferenza(0,0,RSI,m2,q2);

slot_open_ang=acs*2*pi/(ns*p)/2;
[x1,y1]=intersezione_retta_circonferenza(0,0,RSI,tan(slot_open_ang),0);
if (tan(y1./x1)>tan(y1t./x1t))
    x1=x1t;
    y1=y1t;
end

x2=x1+ttd;
y2=y1;

mtta=tan(pi/2-tta*pi/180);
qtta=y2-mtta*x2;
% ytta=mtta*xx+qtta;
[x3,y3]=intersezione_tra_rette(mtta,-1,qtta,m2,-1,q2);

% end of the slot
x6=RSI+lt;
y6=0;
% LT2 position at the tooth
% q2 = y3;
% m2 = 0;
if geo.parallel_slot
    q2 = y3;
    m2 = 0;
end
[xLT2,yLT2]=intersezione_retta_circonferenza(0,0,(RSI+lt),m2,q2);
% bottom slot radius
[xRacSlot,yRacSlot,x5,y5,x4,y4]=cir_tg_2rette(x6,y6,xLT2,yLT2,xLT2,yLT2,x3,y3,RaccordoFC);
x4=x4(2);
y4=y4(2);
x5=x5(2);
y5=y5(2);
xRacSlot=xRacSlot(2);
yRacSlot=yRacSlot(2);

%% slot angle for SPM subdomain model
geo.SlotAngle = atan(y4/x4)*2;                                              % bsa
geo.SlotOpenAngle = slot_open_ang*2;                                        % boa

mm1 = (yLT2-y6)./(xLT2-x6);       % slope of line slot bottom
mm2 = (y3-yLT2)./(x3-xLT2);       % slope of line slot side
angle1 = atan(abs((mm2-mm1)./(1+mm1.*mm2)));        % angle between two lines
%% Chao limit fillet radius 
if x5 > x6    % not beyond central of slot at slot bottom
    x5 = x6;
    y5 = y6;
    RaccordoFC = tan(angle1/2)*sqrt((yLT2-y5)^2+(xLT2-x5)^2);
    [xRacSlot,yRacSlot,x5,y5,x4,y4]=cir_tg_2rette(x6,y6,xLT2,yLT2,xLT2,yLT2,x3,y3,RaccordoFC);
    x4=x4(2);
    y4=y4(2);
    x5=x5(2);
    y5=y5(2);
    xRacSlot=xRacSlot(2);
    yRacSlot=yRacSlot(2);
end
if x4 < x2+(x6-x2)/2    % not beyond half of slot at slot side
    x4 = x2+(x6-x2)/2;
    y4 = m2*x4+q2;
    RaccordoFC = tan(angle1/2)*sqrt((yLT2-y4)^2+(xLT2-x4)^2);
    [xRacSlot,yRacSlot,x5,y5,x4,y4]=cir_tg_2rette(x6,y6,xLT2,yLT2,xLT2,yLT2,x3,y3,RaccordoFC);
    x4=x4(2);
    y4=y4(2);
    x5=x5(2);
    y5=y5(2);
    xRacSlot=xRacSlot(2);
    yRacSlot=yRacSlot(2);
end
geo.SFR = RaccordoFC;
area_corner = RaccordoFC^2 * (1./tan(angle1/2)-(pi-angle1)/2);           % redundant area at bottom slot area (made by slot side, slot bottom and circle)
% slot air (near air-gap)
xA1=RSI;
yA1=0;
% slot air near copper
xA2=x2; yA2=0;
% external point
xE1=RSE; yE1=0;
xE2=RSE*cos(alpha_slot/2); yE2=RSE*sin(alpha_slot/2);
%% slot area evaluation
xArea=[xA2,x2,x3,xLT2,x6,xA2];
yArea=[yA2,y2,y3,yLT2,y6,yA2];
geo.Aslot = 2*(polyarea(xArea,yArea)-area_corner);

% Matrix describing lines and arches of half slot
CavaMat=[0 0 x1 y1 x0 y0 1;
    0 0 xA1 yA1 x1 y1 1;
    x1 y1 x2 y2 NaN NaN 0;
    x2 y2 x3 y3 NaN  NaN 0;
    xRacSlot yRacSlot x5 y5 x4 y4  1;
    x5,y5,x6,0,NaN,NaN,0;
    0 0 xE1 yE1 xE2 yE2 1;
    x0,y0,xE2,yE2,NaN,NaN,0;
    ];
if (q<1 && strcmp(geo.slot_layer_pos,'side_by_side'))
    xstr=(x2+x6)/2;
    ystr=y5/2;
    CavaStr=[x3,y3,x4,y4,NaN,NaN,0];
    CavaMat=[CavaMat;CavaStr];
    
    if ((xA2-xA1)<0.5)
        CavaMat=[CavaMat;[x3,y3,x3,0,NaN,NaN,0;x3,0,x6,y6,NaN,NaN,0]];
    else
        CavaMat=[CavaMat;[x2,y2,x2,0,NaN,NaN,0;x2,0,x6,y6,NaN,NaN,0]];
    end
else
    % lines dividing the slot copper in 2 or more layers (bundles)
    Num_stat_strati=size(avv,1);
    CavaStr=[];
    for jk=1:Num_stat_strati
        xstr(jk)=x3+(jk-1)*(x2-x3)+(jk-1)*(x6-x2)/Num_stat_strati;
        ystr(jk)=y3+(y4-y3)/(x4-x3)*(xstr(jk)-x3);
        CavaStr=[CavaStr;[xstr(jk),0,xstr(jk),ystr(jk),NaN,NaN,0]];
    end
    CavaMat=[CavaMat;CavaStr(2:end,:)];
    if ((xA2-xA1)<0.5)
        CavaMat=[CavaMat;[x3,y3,x3,0,NaN,NaN,0]];
    else
        CavaMat=[CavaMat;[x2,y2,x2,0,NaN,NaN,0]];
    end
    
    for jk=1:Num_stat_strati-1
        CavaMat=[CavaMat;[xstr(jk),ystr(jk),xstr(jk+1),ystr(jk+1),NaN,NaN,0]];
    end
    CavaMat=[CavaMat;[xstr(end),ystr(end),x4,y4,NaN,NaN,0]];
end

% MIRROR of the half slot -> entire slot
CavaMatNeg=CavaMat;
CavaMatNeg(:,2)=-CavaMat(:,2);
CavaMatNeg(:,4)=-CavaMatNeg(:,4);
CavaMatNeg(:,6)=-CavaMat(:,6);
CavaMatNeg(:,7)=-CavaMat(:,7);

CavaMat=[CavaMatNeg;CavaMat];

% replicate/rotate one slot to obtain all slots
[~,ncol] = size(CavaMat);
MatRot=[];
for num_cave=1:(Qs-1)
    jk=1;
    CavaRot=[];
    while jk<=(ncol-2)
        [xrot,yrot]=rot_point(CavaMat(:,jk),CavaMat(:,jk+1),num_cave*alpha_slot);
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
x_air_slot=(xA1+x2)/2;
y_air_slot=0;
Air_slot=[x_air_slot,y_air_slot];
Temp=[];
nome_air_slot=[];
for num_cave=1:(Qs-1)
    
    [xtemp,ytemp]=rot_point(Air_slot(:,1),Air_slot(:,2),num_cave*alpha_slot);
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
[xFe_esp_pos,yFe_esp_pos]=rot_point((RSI+RSI+lt)/2,-wt/2,alpha_slot/2);
xFe_esp_neg=xFe_esp_pos;
yFe_esp_neg=-yFe_esp_pos;
Fe_esp=[xFe_esp_pos,yFe_esp_pos;xFe_esp_neg,yFe_esp_neg];
Temp=[];
nome_Fe_esp=[];
for num_cave=1:(Qs-1)
    
    [xtemp,ytemp]=rot_point(Fe_esp(:,1),Fe_esp(:,2),num_cave*alpha_slot);
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
if (q<1 && strcmp(geo.slot_layer_pos,'side_by_side'))
    x_Cu_slot=[xstr,xstr];
    y_Cu_slot=[ystr,-ystr];
else
    x_Cu_piu=[xstr(2:end),RS5];
    x_Cu_min=xstr;
    x_Cu_slot=(x_Cu_piu+x_Cu_min)/2;
    y_Cu_slot=zeros(1,length(ystr));
end
Cu_slot=[x_Cu_slot',y_Cu_slot'];
Temp=[];
nome_Cu_slot=[];

for num_cave=1:(Qs-1)
    [xtemp,ytemp]=rot_point(x_Cu_slot',y_Cu_slot',num_cave*alpha_slot);
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
xFeYoke=(RS5+RSE)/2;
yFeYoke=0;

FeYoke=[xFeYoke,yFeYoke];
Temp=[];
nome_FeYoke=[];

for num_cave=1:(Qs-1)
    
    [xtemp,ytemp]=rot_point(xFeYoke,yFeYoke,num_cave*alpha_slot);
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
if (Qs<geo.p*geo.q*6)
    codBound_periodic = 10;     % 10 odd or even periodicity
else
    codBound_periodic = -10;    % -10 no periodicity, simulate full machine
end

[xBoundRSE1,yBoundRSE1]=rot_point(RSE,0,alpha_slot/4);  % rotazione di 1/4 di passo cava
[xBoundRSE2,yBoundRSE2]=rot_point(RSE,0,-alpha_slot/4); % rotazione di 1/4 di passo cava
BoundRSE=[[xBoundRSE1,yBoundRSE1];[xBoundRSE2,yBoundRSE2]];
Temp=[];
for num_cave=1:Qs-1
    [xtemp,ytemp]=rot_point([xBoundRSE1;xBoundRSE2],[yBoundRSE1;yBoundRSE2],num_cave*alpha_slot);
    Temp=[Temp;xtemp,ytemp];
end

BoundRSE=[BoundRSE;Temp];
clear Temp xtemp ytemp;
BoundRSE=[BoundRSE,codBound_FluxTan*ones(size(BoundRSE,1),1)];

[xBoundLAT1,yBoundLAT1]=rot_point(mean([RSI,RSE]),0,-alpha_slot/2);
[xBoundLAT2,yBoundLAT2]=rot_point(mean([RSI,RSE]),0,(2*Qs-1)*alpha_slot/2);

BoundLAT=[xBoundLAT1,yBoundLAT1,codBound_periodic;xBoundLAT2,yBoundLAT2,codBound_periodic];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output matrixes:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% matrix 1: lines and arches
statore=CavaMat;

% matrix 2: block labels
if acs~=0
    BLKLABELSstat.xy=[Air_slot;Cu_slot;FeYoke];
    BLKLABELSstat.names.air_slot=nome_air_slot';
    BLKLABELSstat.names.Cu_slot=nome_Cu_slot';
    BLKLABELSstat.names.FeYoke=nome_FeYoke';
    BLKLABELSstat.names.legend={'air_slot','Cu_slot','FeYoke'};
else
    BLKLABELSstat.xy=[Cu_slot;FeYoke];
    BLKLABELSstat.names.air_slot=[];
    BLKLABELSstat.names.Cu_slot=nome_Cu_slot';
    BLKLABELSstat.names.FeYoke=nome_FeYoke';
    BLKLABELSstat.names.legend={'air_slot','Cu_slot','FeYoke'};
end
% matrix 3: buondary conditions
BLKLABELSstat.boundary=[BoundLAT;BoundRSE];

% other data
STATOREdati.DatiVari.RSE=RSE;
STATOREdati.DatiVari.RSI=RSI;
STATOREdati.DatiVari.AngPC=alpha_slot;
STATOREdati.DatiVari.Qs=Qs;


