%% 2013/07/02 MG Following script realize the rotor construction:
% determination of arc and line position
% determinazione of central label position
%% rotore and matrix BLKLABELS

switch geo.RotType
%     case 'ISeg_HS'
%         [geo,temp] = nodes_rotor_ISeg_HS(geo);
%         build_matrix_ISeg;
    case 'ISeg'
        [geo,temp] = nodes_rotor_ISeg(geo);
        build_matrix_ISeg;        
    case 'Seg'
        [geo,temp] = nodes_rotor_Seg(geo);
        build_matrix_Seg;
    case 'Circular'
        [geo,temp] = nodes_rotor_Circ(geo);
        build_matrix_Circ;
    case 'Fluid'
        [geo,temp] = nodes_rotor_Fluid(geo);
        build_matrix_Fluid;
end

codBound_FluxTan=0;
codBound_periodic=10;           % Codice relativo alla periodicità
% 0 flux tangenzial
% 10 odd or even periodicity 

r = geo.r;                     % Raggio del rotore al traferro
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
ang_pont0 = geo.ang_pont0;      % Ampiezza dell'angolo (in gradi) da spazzare con  raggio r in modo da ottenre un arco lungo pont0

nmax = geo.nmax;                % Velocità max (rpm) per la valutazione della sollecitazione centrifuga più gravosa (-> ponticelli)

% dimensione limite della mesh
res=fem.res;
res_traf=fem.res_traf;

%%%%%%%%%%%%%%%%%%%
% Reference to BLKLABELS.materials positions:
codMatFe=5;
codMatBar=6;
codMatShaft=7;
codMatAvvRot=8;
BLKLABELS.materials=geo.BLKLABELSmaterials;
%% %%%%%%%%%%%%%
if not(isempty(rotore))
    rotNeg=rotore;
    rotNeg(:,2)=-rotore(:,2);
    rotNeg(:,4)=-rotore(:,4);
    rotNeg(:,6)=-rotore(:,6);
    rotNeg(:,7)=-rotore(:,7);
    rotore=[rotore;rotNeg];
end

if not(isempty(rotore))
    Temp2=[];
    num_poli=1;
    while num_poli<ps
        Temp1=[];
        for kk=1:2:size(rotore,2)-2
            [xtemp,ytemp]=rot_point(rotore(:,kk),rotore(:,kk+1),num_poli*180/p*pi/180);
            Temp1=[Temp1,xtemp,ytemp];
        end
        Temp2=[Temp2;[Temp1,rotore(:,end)]];
        num_poli=num_poli+1;
    end
    rotore=[rotore;Temp2];
    clear Temp1 Temp2 xtemp ytemp;
end
%
% raggio esterno
xre1 = r;        yre1 = 0;
[xre3,yre3] = rot_point(xre1,yre1,(ps-1/2)*180/p*pi/180);
[xre2,yre2] = rot_point(xre1,yre1,-90/p*pi/180);

rotore=[rotore;
    0 0 xre2 yre2 xre3 yre3 1];

% Albero
[xAl1,yAl1] = rot_point(Ar,0,(ps-1/2)*180/p*pi/180);
[xAl2,yAl2] = rot_point(Ar,0,-90/p*pi/180);
rotore=[rotore;
    0,0,xAl2,yAl2,xAl1,yAl1,1];

% Completo il rotore con linee di chiusura laterali
rotore=[rotore;
    0,0,xre2,yre2,NaN,NaN,0;
    0,0,xre3,yre3,NaN,NaN,0];

%
[nrig,ncol] = size(rotore);

%%  %%
%% Determinazione posizioni baricentriche materiali di ROTORE
%%  %%

%% Baricentri barriere di flux
% non ancora implementata direzione di magnetizzazione
xmag=0*ones(length(xc),1); ymag=0*ones(length(xc),1); zmag=0*ones(length(xc),1);
%
if isempty(xc)
    BarCenter=[];
else
BarCenter=[xc',yc',codMatBar*ones(length(xc),1),res*ones(length(xc),1),1*ones(length(xc),1),xmag,ymag,zmag];
end
Temp=[];
kk=1;
if not(isempty(BarCenter))
    while kk<=ps-1
        [xtemp,ytemp]=rot_point(BarCenter(:,1),BarCenter(:,2),kk*180/p*pi/180);
        Temp=[Temp;xtemp,ytemp,BarCenter(:,3:end)];
        kk=kk+1;
    end
    BarCenter=[BarCenter;Temp];
    clear Temp xtemp ytemp;
end
% Addition negative Barcenter if rad ribs are present
xc_2=[];yc_2=[];
for kk=1:length(xc)
   if  (YpontRadSx(kk)~=0)
      xc_2=[xc_2;xc(kk)]; yc_2=[yc_2;-yc(kk)];
   end
end

if (sum(YpontRadSx~=0)~=0)
    Temp=[];
    kk=1;
    for kk=0:ps-1
        [xtemp,ytemp]=rot_point(xc_2,yc_2,kk*180/p*pi/180);
        Temp=[Temp;xtemp,ytemp,codMatBar*ones(size(xc_2,1),1),res*ones(size(xc_2,1),1),1*ones(size(xc_2,1),1),zeros(size(xc_2,1),1),zeros(size(xc_2,1),1),zeros(size(xc_2,1),1)];
        kk=kk+1;
    end
    
    BarCenter=[BarCenter;Temp];
    clear Temp xtemp ytemp;
end
% Assigment of label name
YpontRadSx=repmat(YpontRadSx,1,nlay*ps);
BarName=[];
kkk=1;
if isempty(BarCenter)
    BarName=[];
else
    for kk=1:length(xc)*ps
        
        if (YpontRadSx(kk)~=0)
            num_bar_element=2;
        else
            num_bar_element=1;
        end
        
        for jk=1:num_bar_element
            if (kk<=nlay)
                ii=kk;
                ij=jk;
            else
                ii=kk-nlay;
                ij=num_bar_element+jk;
            end
            BarName{kkk}={['Air_Bar_',num2str(ii),'_',num2str(ij)]};
            kkk=kkk+1;
        end
    end
end
%% Baricentro ferro di rotore
RotBaricentro=[mean([Ar,XBan1sx(nlay)]),0,codMatFe,res,1,NaN,NaN,NaN];
RotName={'rotor'};

%% Baricentro albero
ShaftBaricentro=[mean([0,Ar]),0,codMatShaft,res,1,NaN,NaN,NaN];
ShaftName={'shaft'};

%% ASSEGNAZIONE BOUNDARY CONDITION
% shaft boundary
[xShaftBound1,yShaftBound1]=rot_point(mean([0,Ar]),0,-90/p*pi/180);
[xShaftBound2,yShaftBound2]=rot_point(mean([0,Ar]),0,(ps-1/2)*180/p*pi/180);
% rotor boundary
[xRotBound1,yRotBound1]=rot_point(mean([Ar,r]),0,-90/p*pi/180);
[xRotBound2,yRotBound2]=rot_point(mean([Ar,r]),0,(ps-1/2)*180/p*pi/180);

clear Mat
%% Matrice centri
BLKLABELSrot.xy=[BarCenter;RotBaricentro;ShaftBaricentro];
BLKLABELSrot.names.BarName=BarName';
BLKLABELSrot.names.RotName=RotName;
BLKLABELSrot.names.ShaftName=ShaftName;
BLKLABELSrot.names.legend={'BarName','RotName','ShaftName'};

BLKLABELSrot.boundary=[xShaftBound1,yShaftBound1,codBound_periodic;
    xShaftBound2,yShaftBound2,codBound_periodic;
    xRotBound1,yRotBound1,codBound_periodic;
    xRotBound2,yRotBound2,codBound_periodic];
%%  %%
%% Rotazione per convenzione segni MOGA
%%  %%

rotore2=[];
for ii=1:2:ncol-2
    [xtemp,ytemp]=rot_point(rotore(:,ii),rotore(:,ii+1),90/p*pi/180);
    rotore2=[rotore2,xtemp,ytemp];
end
rotore2=[rotore2,rotore(:,ncol)];

%% Ruota matrice centri per configurazione rotore MOGA
[xtemp,ytemp]=rot_point(BLKLABELSrot.xy(:,1),BLKLABELSrot.xy(:,2),90/p*pi/180);
BLKLABELSrot.xy=[xtemp,ytemp,BLKLABELSrot.xy(:,3:end)];
clear xtemp ytemp;

%% Rotazione boundary xy
[xtemp,ytemp]=rot_point(BLKLABELSrot.boundary(:,1),BLKLABELSrot.boundary(:,2),90/p*pi/180);
BLKLABELSrot.boundary=[xtemp,ytemp,BLKLABELSrot.boundary(:,3:end)];
%%
