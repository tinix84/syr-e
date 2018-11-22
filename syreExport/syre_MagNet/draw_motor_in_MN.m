

function [geo,mat] = draw_motor_in_MN(geo,eval_type,mat,pathname,filename,h)

%% tmp folder for MagNet export
% pathname=[cd '\tmp\MagnetExportTmp\']; % creo cartella temporanea
% mkdir(pathname);
% filename='mot_0.mat';

fem = dimMesh(geo,eval_type);

% calc winding factor (kavv) and rotor offset (phase1_offset)
[kavv, phase1_offset] = calcKwTh0(geo.tempWinTable,geo.ns*geo.p,geo.p);

% offset angle for coordinate transformations
if strcmp(geo.RotType,'SPM')
    phase1_offset = phase1_offset-90;
end
th_m0 = 0;                              % rotor position [mec deg]
geo.th0 = th_m0*geo.p - phase1_offset;  % d- to alpha-axis offset [elt deg]

% nodes
geo.x0 = geo.r/cos(pi/2/geo.p);
geo.fem=fem;

[rotor,BLKLABELSrot,geo,mat] = ROTmatr(geo,fem,mat); % rotor and BLKLABELSrot describe the rotor
geo.rotor = rotor;

[geo,statore,BLKLABELSstat] = STATmatr(geo,fem); % statore and BLKLABELSstat describe the stator
geo.stator=statore;

% CodificaMateriali;
% CodificaMateriali_prova;
% BLKLABELS.materials=char(Air1, Air2, Stat_Cu,Stat_Fe,Rot_Fe,Stat_Magn,Rot_Magn,Shaft_Steel,Rot_Cu);  %scelgo io i nomi dei materiali
% BLKLABELS.materials=char(geo.BLKLABELSmaterials);

BLKLABELS.materials=char('AIR', 'AIR', 'Copper: 5.77e7 Siemens/meter','M600-50A','M600-50A','PM04: Br 0.4 mur 1.0','CR10: Cold rolled 1010 steel','Copper: 5.77e7 Siemens/meter');  %scelgo io i nomi dei materiali
% BLKLABELS.materials=char('AIR', 'AIR', 'Copper: 5.77e7 Siemens/meter','M250-35A','M250-35A','FB14H','CR10: Cold rolled 1010 steel','Copper: 5.77e7 Siemens/meter');  %scelgo io i nomi dei materiali

BLKLABELS.rotore = BLKLABELSrot;
BLKLABELS.statore= BLKLABELSstat;
geo.BLKLABELS=BLKLABELS;

syreToDxf(geo.stator, geo.rotor,pathname, filename);

% inizio funzione vecchia
CGD='Call getDocument()';
CGDGV='Call getDocument().getView()';
IMCUS='infoMakeComponentUnionSurfaces';
IMCRV='infoMakeComponentRemoveVertices';
ITIS=' infoToggleInSelection';
IATS='infoAddToSelection';
% if ~strcmp(geo.RotType,'SPM')
%     if geo.BarFillFac~=0
%         BarNum=length(BLKLABELSrot.BarName);
%         for ii=1:BarNum/2
%             BLKLABELSrot.BarName{ii+BarNum}={'Barrier_Air_' int2str(ii)};
%         end
%     end
% end

Br = mat.LayerMag.Br(1);
% variable not present in geo
Q=geo.ns*geo.p;      %Mac.Q=6*Mac.p*Mac.q;
%Mac.RtS=geo.r+geo.g; (non usata)
%Mac.Br=mat.LayerMag.Br; (non usata)
% ========================
l=geo.l;

% h = OpenMagnet(1);  % 1 significa visibile, 0 invisibile

% Nuovo Documento
h = NewDocumentMagnet(h);

% h = SaveDocumentMagnet(h,docName);

% Importa file dxf e crea geometria
% ImportDXFMagnet(h,dxfName);
ImportDXFMagnet(h,[pathname filename(1:end-4),'.dxf']);

%% prova confronto materiali
% mh = h.magnetHandler;
% Command = ['CALL getDocument().getModelMaterialDatabase().getMaterials(materials)']
% invoke(mh, 'processCommand', Command);

% MN6 = actxserver('Magnet.application');
% Command = ['CALL getSystemMaterialDatabase().getMaterials(BMN-38,EH)']
% invoke(MN6, 'processCommand', Command);
%%

pc=360/(Q)/2; % mezzo passo cava

[traferro,BLKLABELStraf]=traferroMatr(pc,geo.r,geo.g,Q,geo.Qs,geo.ps,geo.p,fem.res_traf);
[rig_traf, col_traf] = size(traferro);
if strcmp(geo.RotType,'SPM')
    % SF (30/08/2018): for the SPM machines the airgap cover also the air
    % space on the PMs side, so the lines must be extended
    tmp=traferro(5,:);
    tmpAng=atan2(tmp(2),tmp(1));
    tmpRad=geo.r-geo.lm;
    tmp(1)=tmpRad*cos(tmpAng);
    tmp(2)=tmpRad*sin(tmpAng);
    traferro(5,:)=tmp;
    tmp=traferro(6,:);
    tmpAng=atan2(tmp(2),tmp(1));
    tmpRad=geo.r-geo.lm;
    tmp(1)=tmpRad*cos(tmpAng);
    tmp(2)=tmpRad*sin(tmpAng);
    traferro(6,:)=tmp;
end

for ii=1:rig_traf
    if(traferro(ii,col_traf)==0)
        DrawLineMagnet(h,[traferro(ii,1) traferro(ii,2)],[traferro(ii,3) traferro(ii,4)]);
    else
        DrawArcMagnetperPunti(h,traferro(ii,1),traferro(ii,2),traferro(ii,3),traferro(ii,4),traferro(ii,5),traferro(ii,6));
    end
end

% Assegna materiali alle diverse parti della macchina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statore
[var1,var2]=size(BLKLABELSstat.names.air_slot);
for kk=1:var1
    h = MakeComponentMagnet(h,[BLKLABELSstat.xy(kk,1), BLKLABELSstat.xy(kk,2)],cell2mat(BLKLABELSstat.names.air_slot{kk,1}),l,BLKLABELS.materials(2,:),'None',[0, 0, 1],BLKLABELSstat.xy(kk,4));
    indice_air_slot=kk;
end
[var1,var2]=size(BLKLABELSstat.names.Cu_slot);
for kk=1:var1
    h = MakeComponentMagnet(h,[BLKLABELSstat.xy(kk+indice_air_slot,1), BLKLABELSstat.xy(kk+indice_air_slot,2)],cell2mat(BLKLABELSstat.names.Cu_slot{kk,1}),l,BLKLABELS.materials(3,:),'None',[0, 0, 1],-1);
    indice_Cu_slot=kk+indice_air_slot;
end
[var1,var2]=size(BLKLABELSstat.names.FeYoke);
for kk=1:var1
    h = MakeComponentMagnet(h,[BLKLABELSstat.xy(kk+indice_Cu_slot,1), BLKLABELSstat.xy(kk+indice_Cu_slot,2)],cell2mat(BLKLABELSstat.names.FeYoke{kk,1}),l,BLKLABELS.materials(4,:),'None',[0, 0, 1],BLKLABELSstat.xy(kk+indice_Cu_slot,4));
end

% Rotore
for kk=1:length(BLKLABELSrot.xy(:,1))
    switch BLKLABELSrot.xy(kk,3)
        case 1 %air
            h = MakeComponentMagnet(h,[BLKLABELSrot.xy(kk,1), BLKLABELSrot.xy(kk,2)],(BLKLABELSrot.BarName{kk,1}),l,BLKLABELS.materials(2,:),'Uniform',[BLKLABELSrot.xy(kk,6),BLKLABELSrot.xy(kk,7), BLKLABELSrot.xy(kk,8)],BLKLABELSrot.xy(kk,4));
        case 5 %rotor iron
            h = MakeComponentMagnet(h,[BLKLABELSrot.xy(kk,1), BLKLABELSrot.xy(kk,2)],'rotor',l,BLKLABELS.materials(5,:),'None',[0, 0, 1],BLKLABELSrot.xy(kk,4));
        case 6 %PM
            if Br==0 % || strcmp(geo.RotType,'SPM')
                h = MakeComponentMagnet(h,[BLKLABELSrot.xy(kk,1), BLKLABELSrot.xy(kk,2)],(BLKLABELSrot.BarName{kk,1}),l,BLKLABELS.materials(2,:),'None',[0, 0, 1],BLKLABELSrot.xy(kk,4));
            else
                %MatName=BLKLABELS.materials(6,:);
                h = MakeComponentMagnet(h,[BLKLABELSrot.xy(kk,1), BLKLABELSrot.xy(kk,2)],(BLKLABELSrot.BarName{kk,1}),l,BLKLABELS.materials(6,:),'Uniform',[BLKLABELSrot.xy(kk,6),BLKLABELSrot.xy(kk,7), BLKLABELSrot.xy(kk,8)],BLKLABELSrot.xy(kk,4));
            end
        case 7 %shaft
            h = MakeComponentMagnet(h,[BLKLABELSrot.xy(kk,1), BLKLABELSrot.xy(kk,2)],'shaft',l,BLKLABELS.materials(7,:),'None',[0, 0, 1],-1);
    end
end
% Traferro
for kk=1:size(BLKLABELStraf.names,1)
    h = MakeComponentMagnet(h,[BLKLABELStraf.xy(kk,1), BLKLABELStraf.xy(kk,2)],BLKLABELStraf.names{kk,1},l,BLKLABELS.materials(BLKLABELStraf.xy(kk,3),:),'None',[0, 0, 1],-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Crea le bobine di statore
avv=geo.avv(:,1:geo.Qs);

% COIL U
coil_name = 'U';
coil_number = 1;
MakeSimpleCoilMagnet_UVW;

%COSTRUZIONE COIL V
coil_name = 'V';
coil_number = 2;
MakeSimpleCoilMagnet_UVW;

%COSTRUZIONE COIL W
coil_name = 'W';
coil_number = 3;
MakeSimpleCoilMagnet_UVW

% Cancella tutte le linee di costruzione
DeleteAllLinesMagnet(h, geo.R+50, -geo.R-50, -geo.R-50, geo.R+50);

%%%%%%%%%%%%%%%%%%
% Costruzione delle airbox
%%%%%%%%%%%%%%%%%%
%% 2013/07/29 MG in questa prima bozza si esegue una rivalutazione dei punti principali al traferro
[xArcTraf1,yArcTraf1] = rot_point(geo.r+1/3*geo.g,0,-pc*pi/180);
[xArcTraf2,yArcTraf2] = rot_point(geo.r+1/3*geo.g,0,(2*geo.Qs-1)*pc*pi/180);

[THETA1,RHO1] = cart2pol(xArcTraf1,yArcTraf1);
THETA1=THETA1*180/pi;
[THETA2,RHO2] = cart2pol(xArcTraf2,yArcTraf2);  %RHO2 non usato
THETA2=THETA2*180/pi;

DrawArcSectorMagnet(h,[0 0],RHO1,geo.R,THETA1,THETA2);
DrawArcMagnet(h,[0 0],RHO1,THETA2,geo.ps*180/geo.p);

DrawLineMagnet(h,[0 0],[RHO1*cos(geo.ps*pi/geo.p) RHO1*sin(geo.ps*pi/geo.p)]);
DrawLineMagnet(h,[0 0],[RHO1 0]);

h = MakeComponentMagnet(h,[0, 0],'Rotor_Airbox',l,BLKLABELS.materials(1,:),'None',[0, 0, 1],-1);
h = MakeComponentMagnet(h,[RHO1*1.001, RHO1*0.0001],'Stator_Airbox',l,BLKLABELS.materials(1,:),'None',[0, 0, 1],-1);

% Suddivisione mesh al traferro
clear edgeList;
edgeList{1}{1} = 'air_gap_rotor';
edgeList{1}{2} = 2;
edgeList{2}{1} = 'air_gap_stator1';
edgeList{2}{2} = 2;
edgeList{3}{1} = 'air_gap_stator2';
edgeList{3}{2} = 2;

% (MG) Machine periodicity selection
% t = gcd(round(Q),geo.p);  % number of periods
if rem(geo.ps,2)==0
    bdryType='Even';
else
    bdryType='Odd';
end

if (geo.ps*180/geo.p<=90)
    edgeList{1}{3} = 3;
    edgeList{2}{3} = 3;
    edgeList{3}{3} = 3;
else
    edgeList{1}{3} = 4;
    edgeList{2}{3} = 4;
    edgeList{3}{3} = 4;
end

DimMesh=round(2*pc*geo.Qs/fem.res_traf);
SetUniformMeshOnEdgesMagnet(h,edgeList,DimMesh);

%================ Imposta condizioni al contorno ==========================

% Flux Tangential su superficie esterna di Statore
clear bdryFaces
bdryFaces{1}{1} = 'Stator_Airbox';
if (geo.ps*180/geo.p<=90)
    bdryFaces{1}{2} = 5;
else
    bdryFaces{1}{2} = 6;
end
bdryName = 'Flux Tangential';
SetBdryTangentialMagnet(h,bdryName,bdryFaces);

% Odd or even boundary aplication
clear bdryFaces
if (geo.ps*180/geo.p<=90)
    bdryFaces{1}{1} = 'Stator_Airbox';
    bdryFaces{1}{2} = 6;
    bdryFaces{2}{1} = 'Rotor_Airbox';
    bdryFaces{2}{2} = 5;
    bdryFaces{3}{1} = 'Rotor_Airbox';
    bdryFaces{3}{2} = 4;
elseif (geo.ps*180/geo.p==180)
    bdryFaces{1}{1} = 'Stator_Airbox';
    bdryFaces{1}{2} = 3;
    bdryFaces{2}{1} = 'Rotor_Airbox';
    bdryFaces{2}{2} = 4;
    bdryFaces{3}{1} = 'Rotor_Airbox';
    bdryFaces{3}{2} = 3;
else
    bdryFaces{1}{1} = 'Stator_Airbox';
    bdryFaces{1}{2} = 3;
    bdryFaces{2}{1} = 'Rotor_Airbox';
    bdryFaces{2}{2} = 5;
    bdryFaces{3}{1} = 'Rotor_Airbox';
    bdryFaces{3}{2} = 3;
end

bdryName = 'Periodic';
rotAxis = [0, 0, 1];
rotAngle = -(geo.ps*180/geo.p);
SetBdryRoundPeriodicMagnet(h,bdryName,bdryType,bdryFaces,rotAxis,[0, 0, 0],rotAngle)

% Imposta parte rotante

[var1 var2]=size(BLKLABELSrot.BarName);     %var2
clear Motion_Components
for kk=1:var1
    Motion_Components{kk}=(BLKLABELSrot.BarName{kk,1});
end
if strcmp(geo.RotType,'SPM')
    %     for jj = kk+1:kk+geo.ps*2
    %         Motion_Components{jj}=['Rotor_Air_Zone_',num2str(jj-kk)];
    %     end
    
    jj = 1;
    Motion_Components{jj+1}='rotor';
    Motion_Components{jj+2}='shaft';
    Motion_Components{jj+3}= 'air_gap_rotor';
    Motion_Components{jj+4}= 'Rotor_Airbox';
else
    Motion_Components{kk+1}='rotor';
    Motion_Components{kk+2}='shaft';
    Motion_Components{kk+3}= 'air_gap_rotor';
    Motion_Components{kk+4}= 'Rotor_Airbox';
end
h = MakeMotionComponentMagnet(h,Motion_Components,'Motion#1',6,0);

%% Imposta Parametri per simulazione Transient with Motion

% SetTransientOptions (h,num2str(tStart),num2str(tStep),num2str(tStop),'Yes');
SetTransientOptions (h,0,0.1,6,'Yes',0)

%% ===================== Salva e chiudi MagNet ============================
% pathname_DXF=[pathname,filename(1:end-4),'\DXF\'];
% pathname_DXF=[pathname,filename(1:end-4)];

% [h,f] = SaveDocumentMagnet(h,[pathname,filename(1:end-4),'.mn']);
% CloseMagnet(h)


