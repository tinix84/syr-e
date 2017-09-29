% Copyright 2015
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

% BuildMachine.m - builds the machine model into Magnet by Infolytica
% input:  motorname.mat (made by syre) and motorname.dxf (from syreToDxf.m)
% output: motorname.mn in the same folder

clear all, close all, clc
%%%%%%%%%%%%%%%%%%%

p = genpath('MNscripts_dxfBuild')
addpath(p);

%% Commonly used words
CGD='Call getDocument()';
CGDGV='Call getDocument().getView()';
IMCUS='infoMakeComponentUnionSurfaces';
IMCRV='infoMakeComponentRemoveVertices';
ITIS=' infoToggleInSelection';
IATS='infoAddToSelection';

%% === Apertura di Magnet e creazione modello =============================
load ultimo;
[docName,pathname]=uigetfile([pathname,'\*.mat'],'Select the Machine Mat-file');
save ultimo.mat pathname -append;
load([pathname,docName]);
dxfName=[pathname,docName(1:end-4),'.dxf'];
% CENTRI=BLKLABELS;

CENTRIstat=geo.BLKLABELS.statore;
CENTRIrot=geo.BLKLABELS.rotore;

if ~strcmp(geo.RotType,'SPM')
    if geo.BarFillFac~=0
        BarNum=length(CENTRIrot.BarName);
        for ii=1:BarNum/2
            CENTRIrot.BarName{ii+BarNum}={'Barrier_Air_' int2str(ii)};
        end
    end
end
CodificaMateriali;
CENTRI.materials=char(Air1, Air2, Stat_Cu,Stat_Fe,Rot_Fe,Stat_Magn,Rot_Magn,Shaft_Steel,Rot_Cu);
% CENTRI.materials=char(Air1, Air2, Stat_Cu,Stat_Fe,Rot_Fe,Rot_Magn,Shaft_Steel,Rot_Cu);
Mac=geo;
Mat = mat;
Br = mat.LayerMag.Br;
%% variable not present in geo
Mac.Q=6*Mac.p*Mac.q;
Mac.RtS=geo.r+geo.g;
Mac.Br=mat.LayerMag.Br;
%% ========================
L_assiale=Mac.l;
t=gcd(Mac.Q,Mac.p); % machine periodicity
% keyboard
h = OpenMagnet(1);  % 1 significa visibile, 0 invisibile

% Nuovo Documento
h = NewDocumentMagnet(h);

% h = SaveDocumentMagnet(h,docName);

% Importa file dxf e crea geometria
ImportDXFMagnet(h,dxfName);

pc=360/(Mac.Q)/2; % mezzo passo cava
xr=Mac.RtS-Mac.g;

[traferro,CENTRItraf]=traferroMatr(pc,xr,Mac.g,Mac.Q,Mac.Qs,Mac.ps,Mac.p,Mac.fem.res_traf);
[rig_traf, col_traf] = size(traferro);

% Se geometria MOGA
for ii=1:rig_traf
    if(traferro(ii,col_traf)==0)
        DrawLineMagnet(h,[traferro(ii,1) traferro(ii,2)],[traferro(ii,3) traferro(ii,4)]);
        %         keyboard
    else
        DrawArcMagnetperPunti(h,traferro(ii,1),traferro(ii,2),traferro(ii,3),traferro(ii,4),traferro(ii,5),traferro(ii,6));
        %         keyboard
    end
end
% % Se geometria MOGA
% DrawLineMagnet(h,[r_rot*cosd(360/Ksym) r_rot*sind(360/Ksym)],[r1*cosd(360/Ksym) r1*sind(360/Ksym)]);
% DrawLineMagnet(h,[r_rot 0],[r1 0]);
% DrawArcMagnet(h,[0 0],r1,0,360/Ksym);
% DrawArcSectorMagnet(h,[0 0],r2,r3,0-th_Stat_Rot,360/Ksym-th_Stat_Rot);
% DrawArcSectorMagnet(h,[0 0],r1,r2,0-th_Stat_Rot,360/Ksym-th_Stat_Rot);
% % DrawArcMagnet(h,[0 0],r3,0-th_Stat_Rot,360/Ksym-th_Stat_Rot);
% % DrawArcMagnet(h,[0 0],r2,0-th_Stat_Rot,360/Ksym-th_Stat_Rot);

%% Crea componenti nel modello

% Assegna materiali alle diverse parti della macchina

%%%%%%%%%%%%%
%% Statore
%%%%%%%%%%%%%

[var1,var2]=size(CENTRIstat.names.air_slot);
for kk=1:var1
    h = MakeComponentMagnet(h,[CENTRIstat.xy(kk,1), CENTRIstat.xy(kk,2)],cell2mat(CENTRIstat.names.air_slot{kk,1}),L_assiale,CENTRI.materials(2,:),'None',[0, 0, 1],CENTRIstat.xy(kk,4));
    indice_air_slot=kk;
end

[var1,var2]=size(CENTRIstat.names.Cu_slot);
for kk=1:var1
    h = MakeComponentMagnet(h,[CENTRIstat.xy(kk+indice_air_slot,1), CENTRIstat.xy(kk+indice_air_slot,2)],cell2mat(CENTRIstat.names.Cu_slot{kk,1}),L_assiale,CENTRI.materials(3,:),'None',[0, 0, 1],-1);
    indice_Cu_slot=kk+indice_air_slot;
end
[var1,var2]=size(CENTRIstat.names.FeYoke);
for kk=1:var1
    h = MakeComponentMagnet(h,[CENTRIstat.xy(kk+indice_Cu_slot,1), CENTRIstat.xy(kk+indice_Cu_slot,2)],cell2mat(CENTRIstat.names.FeYoke{kk,1}),L_assiale,CENTRI.materials(4,:),'None',[0, 0, 1],CENTRIstat.xy(kk+indice_Cu_slot,4));
end

%%%%%%%%%%%
%% Rotore
%%%%%%%%%%%
% Barrier
% Mac.Br=0.4;
% CENTRI.materials(7,:)='NMF3F';
if Br==0 || strcmp(geo.RotType,'SPM')
    [BarSize,var2]=size(CENTRIrot.BarName);
    for kk=1:BarSize
        h = MakeComponentMagnet(h,[CENTRIrot.xy(kk,1), CENTRIrot.xy(kk,2)],cell2mat(CENTRIrot.BarName{kk,1}),L_assiale,CENTRI.materials(2,:),'None',[0, 0, 1],CENTRIrot.xy(kk,4));
    end
else
    %Barrier Filled with magnet
    [BarSize,var2]=size(CENTRIrot.BarName);
    [LabelSize,~]=size(CENTRIrot.xy);
    BarSize=LabelSize-2;
    for kk=1:LabelSize-2
        if CENTRIrot.xy(kk,3)==6
            tmpmat=CENTRI.materials(7,:);
        else
            tmpmat='AIR';
        end
        h = MakeComponentMagnet(h,[CENTRIrot.xy(kk,1), CENTRIrot.xy(kk,2)],cell2mat(CENTRIrot.BarName{kk,1}),L_assiale,tmpmat,'Uniform',[CENTRIrot.xy(kk,6),CENTRIrot.xy(kk,7), CENTRIrot.xy(kk,8)],CENTRIrot.xy(kk,4));
    end
    
    % for ii=1:length(CENTRIrot.xy(:,1))
    %     switch CENTRIrot.xy(ii,3)
    %         case 1 %air
    %             MatName=CENTRI.materials(1,:);
    %             PartName='AIR';
    %             magType='None';
    %         case 5 %rotor iron
    %             MatName=CENTRI.materials(5,:);
    %             PartName='rotor';
    %             magType='None';
    %             CENTRIrot.xy(6:8)=[0 0 1]; % to avoid NaN in MagNet
    %         case 6 %PM
    %             MatName=CENTRI.materials(7,:);
    %             PartName=cell2mat(CENTRIrot.BarName{kk,1});
    %             magType='Uniform';
    %         case 7 %shaft
    %             MatName=CENTRI.materials(8,:);
    %             PartName='shaft';
    %             magType='None';
    %             CENTRIrot.xy(6:8)=[0 0 1]; % to avoid NaN in MagNet
    %     end
    %
    %     h = MakeComponentMagnet(h,[CENTRIrot.xy(kk,1), CENTRIrot.xy(kk,2)],PartName,L_assiale,MatName,magType,[CENTRIrot.xy(kk,6),CENTRIrot.xy(kk,7), CENTRIrot.xy(kk,8)],CENTRIrot.xy(kk,4));
    %
    % end
end

if strcmp(geo.RotType,'SPM')
%     for jj = kk+1:kk+geo.ps*2
%         h = MakeComponentMagnet(h,[CENTRIrot.xy(jj,1), CENTRIrot.xy(jj,2)],['Rotor_Air_Zone_',num2str(jj-kk)],L_assiale,CENTRI.materials(1,:),'None',[0, 0, 1],-1);
%     end
    
    jj = 1;
    % rotor Iron
    jj=jj+1;
    h = MakeComponentMagnet(h,[CENTRIrot.xy(jj,1), CENTRIrot.xy(jj,2)],'rotor',L_assiale,CENTRI.materials(5,:),'None',[0, 0, 1],CENTRIrot.xy(jj,4));
    
    % Shaft
    jj=jj+1;
    h = MakeComponentMagnet(h,[CENTRIrot.xy(jj,1), CENTRIrot.xy(jj,2)],'shaft',L_assiale,CENTRI.materials(8,:),'None',[0, 0, 1],-1);
else
    % rotor Iron
    kk= kk+1;
    h = MakeComponentMagnet(h,[CENTRIrot.xy(kk,1), CENTRIrot.xy(kk,2)],'rotor',L_assiale,CENTRI.materials(5,:),'None',[0, 0, 1],CENTRIrot.xy(kk,4));
    
    % Shaft
    kk= kk+1;
    h = MakeComponentMagnet(h,[CENTRIrot.xy(kk,1), CENTRIrot.xy(kk,2)],'shaft',L_assiale,CENTRI.materials(8,:),'None',[0, 0, 1],-1);
end

%%%%%%%%%%%%
%% Traferro
%%%%%%%%%%%%
for kk=1:size(CENTRItraf.names,1)
    h = MakeComponentMagnet(h,[CENTRItraf.xy(kk,1), CENTRItraf.xy(kk,2)],CENTRItraf.names{kk,1},L_assiale,CENTRI.materials(CENTRItraf.xy(kk,3),:),'None',[0, 0, 1],-1);
end


%%
% % Crea le bobine di statore
% stringhe di uso comune

avv=Mac.avv(:,1:Mac.Qs);

%COSTRUZIONE COIL U
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

%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%
%% Cancella tutte le linee di costruzione

DeleteAllLinesMagnet(h, Mac.R+50, -Mac.R-50, -Mac.R-50, Mac.R+50);

%%%%%%%%%%%%%%%%%%
% Costruzione delle airbox
%%%%%%%%%%%%%%%%%%
%% 2013/07/29 MG in questa prima bozza si esegue una rivalutazione dei punti principali al traferro
[xArcTraf1,yArcTraf1] = rot_point(xr+1/3*Mac.g,0,-pc*pi/180);
[xArcTraf2,yArcTraf2] = rot_point(xr+1/3*Mac.g,0,(2*Mac.Qs-1)*pc*pi/180);

[THETA1,RHO1] = cart2pol(xArcTraf1,yArcTraf1);
THETA1=THETA1*180/pi;
[THETA2,RHO2] = cart2pol(xArcTraf2,yArcTraf2);
THETA2=THETA2*180/pi;

DrawArcSectorMagnet(h,[0 0],RHO1,Mac.R,THETA1,THETA2);
DrawArcMagnet(h,[0 0],RHO1,THETA2,Mac.ps*180/Mac.p);

DrawLineMagnet(h,[0 0],[RHO1*cos(Mac.ps*pi/Mac.p) RHO1*sin(Mac.ps*pi/Mac.p)]);
DrawLineMagnet(h,[0 0],[RHO1 0]);

h = MakeComponentMagnet(h,[0, 0],'Rotor_Airbox',L_assiale,CENTRI.materials(1,:),'None',[0, 0, 1],-1);
h = MakeComponentMagnet(h,[RHO1*1.001, RHO1*0.0001],'Stator_Airbox',L_assiale,CENTRI.materials(1,:),'None',[0, 0, 1],-1);

% Suddivisione mesh al traferro
clear edgeList;
edgeList{1}{1} = 'air_gap_rotor';
edgeList{1}{2} = 2;
edgeList{2}{1} = 'air_gap_stator1';
edgeList{2}{2} = 2;
edgeList{3}{1} = 'air_gap_stator2';
edgeList{3}{2} = 2;

% (MG) Machine periodicity selection
Q = Mac.Q;                    % number of slots
t = gcd(round(Mac.Q),Mac.p);  % number of periods
% if ((6*t/Q)>1)
%     bdryType='Even';
% else
%     bdryType = 'Odd';
% end

% (SF) Machine periodicity
if rem(Mac.ps,2)==0
    bdryType='Even';
else
    bdryType='Odd';
end

if (Mac.ps*180/Mac.p<=90)
    edgeList{1}{3} = 3;
    edgeList{2}{3} = 3;
    edgeList{3}{3} = 3;
else
    edgeList{1}{3} = 4;
    edgeList{2}{3} = 4;
    edgeList{3}{3} = 4;
end

DimMesh=round(2*pc*Mac.Qs/Mac.fem.res_traf);
SetUniformMeshOnEdgesMagnet(h,edgeList,DimMesh);

%================ Imposta condizioni al contorno ==========================

% Flux Tangential su superficie esterna di Statore
clear bdryFaces
bdryFaces{1}{1} = 'Stator_Airbox';
if (Mac.ps*180/Mac.p<=90)
    bdryFaces{1}{2} = 5;
else
    bdryFaces{1}{2} = 6;
end
bdryName = 'Flux Tangential';
SetBdryTangentialMagnet(h,bdryName,bdryFaces);

% Odd or even boundary aplication
clear bdryFaces
if (Mac.ps*180/Mac.p<=90)
    bdryFaces{1}{1} = 'Stator_Airbox';
    bdryFaces{1}{2} = 6;
    bdryFaces{2}{1} = 'Rotor_Airbox';
    bdryFaces{2}{2} = 5;
    bdryFaces{3}{1} = 'Rotor_Airbox';
    bdryFaces{3}{2} = 4;
elseif (Mac.ps*180/Mac.p==180)
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
rotAngle = -(Mac.ps*180/Mac.p);
SetBdryRoundPeriodicMagnet(h,bdryName,bdryType,bdryFaces,rotAxis,[0, 0, 0],rotAngle)

% Imposta parte rotante

[var1 var2]=size(CENTRIrot.BarName);
clear Motion_Components
for kk=1:var1
    Motion_Components{kk}=cell2mat(CENTRIrot.BarName{kk,1});
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
[h,f] = SaveDocumentMagnet(h,[pathname,docName(1:end-4),'.mn']);
CloseMagnet(h)

