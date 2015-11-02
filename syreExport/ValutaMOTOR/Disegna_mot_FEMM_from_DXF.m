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

%% 2013/07/05 MG Importa i file dxf di statore e di rotore e assegna
%% materiali e condizioni al contorno
% Input: viene caricato il file dxf di cui di vuole ricavare il
% corrispondente file.femm, importante è che il file dxf es motXX.dxf sia
% accompagnato da un file motXX_ParMachine.mat altrimenti il file femm darà
% errore.
% Output: dialogBox richiedente dove scaricare il file femm generato, tale
% file femm sarà accompagnato dal file ParMachine del dxf.

clear all
close all

current_path=cd;
[pathstr, name, ext] = fileparts(current_path);

% addpath C:\Programmi\femm42\mfiles\;
% addpath('C:\Program Files (x86)\femm42\mfiles');
addpath('c:\femm42_beta\mfiles\');
addpath([current_path,'\dxf_conv_fun\']);
addpath([current_path,'\mfiles\']);

error_code = 0;

[filemot, pathname] = uigetfile([cd '\*.dxf'], 'Pick a motor');
load([pathname,filemot(1:end-4),'.mat']);

%% Definizione dei gruppi:
groupStat=1;
groupRot=2;
groupTraf=20;
fem=geo.fem;
%% Definisco la fase meccanica per il momento manualmente e pari a 0 secondo la convenzione di disegno adottata da MOGA:
th_m0=0;

filename = 'mot0.fem';
% opendocument('empty_case.fem',h_temp);
openfemm;
opendocument([current_path,'\empty_case.fem']);
mi_saveas(filename);
mi_probdef(0,'millimeters','planar',1e-8,geo.l,10);

%% Scelta periodicità per i boundary:
if (rem(geo.ps,2)==0)
    Period=4;
else
    Period=5;
end
%% boundary
% A = 0 sulle superfici interna ed esterna
mi_addboundprop('A=0', 0, 0, 0, 0, 0, 0, 0, 0, 0);

% anti-periodicità (2 x rotore + 2 x statore + 3 x airgap + AP move)
mi_addboundprop('APr1', 0, 0, 0, 0, 0, 0, 0, 0, Period);
mi_addboundprop('APr2', 0, 0, 0, 0, 0, 0, 0, 0, Period);

mi_addboundprop('APg1', 0, 0, 0, 0, 0, 0, 0, 0, Period);
mi_addboundprop('APg2', 0, 0, 0, 0, 0, 0, 0, 0, Period);
mi_addboundprop('APg3', 0, 0, 0, 0, 0, 0, 0, 0, Period);
mi_addboundprop('APmove', 0, 0, 0, 0, 0, 0, 0, 0,Period);

mi_addboundprop('APs1', 0, 0, 0, 0, 0, 0, 0, 0, Period);
mi_addboundprop('APs2', 0, 0, 0, 0, 0, 0, 0, 0, Period);

%% DXF file import;
mi_readdxf([pathname,filemot]);
keyboard
BLKLABELS=geo.BLKLABELS;
%% Assign rotor block property 
assign_block_prop_rot(BLKLABELS,geo,fem,2);
% boundary conditions
BLKLABELSrot=BLKLABELS.rotore;
for ii=1:2
    mi_selectsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2));
    if (BLKLABELSrot.boundary(ii,3)==10)
        mi_setsegmentprop('APr1', 0, 1, 0, 2);
        mi_clearselected;
    elseif(BLKLABELSrot.boundary(ii,3)==0)
        mi_selectarcsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2))
        mi_setarcsegmentprop(fem.res, 'A=0', 0, 2);
        mi_clearselected;
    end
end
% Condizioni al contorno di rotore ferro lamierino
for ii=3:4
    
    mi_selectsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2));
    if (BLKLABELSrot.boundary(ii,3)==10)
        mi_setsegmentprop('APr2', 0, 1, 0, 2);
        mi_clearselected;
    elseif(BLKLABELSrot.boundary(ii,3)==0)
        mi_selectarcsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2))
        mi_setarcsegmentprop(fem.res, 'A=0', 0, 2);
        mi_clearselected;
    end
end

%% Assign stator block property
assign_block_prop_stat(BLKLABELS,geo,fem,1) % assegna materiali;
% boundary conditions
BLKLABELSstat=BLKLABELS.statore;
for ii=1:size(BLKLABELSstat.boundary,1)
    
    mi_selectsegment(BLKLABELSstat.boundary(ii,1),BLKLABELSstat.boundary(ii,2));
    if (BLKLABELSstat.boundary(ii,3)==10)
        mi_setsegmentprop('APs1', 0, 1, 0, 1);
        mi_clearselected;
    elseif(BLKLABELSstat.boundary(ii,3)==0)
        mi_selectarcsegment(BLKLABELSstat.boundary(ii,1),BLKLABELSstat.boundary(ii,2));
        mi_setarcsegmentprop(fem.res, 'A=0', 0, 1);
        mi_clearselected;
    end
end

%% Airgap
pc=360/(6*geo.p*geo.q)/2;        % half stator slot
AirGapBuild(geo.Qs,geo.ps,geo.p,geo.g,pc,geo.r,fem.res_traf,1,2);
draw_airgap_arc_with_mesh(geo,th_m0,fem);

%% Save file...
[FileName,PathName,FilterIndex] = uiputfile('*.fem','Save motor in',[pathname,filemot(1:end-4)]);
mi_saveas([PathName,FileName(1:end-4),'.fem']);
closefemm;

% copyfile([pathname,filemot(1:end-4),'_ParMachine.mat'],[PathName,FileName(1:end-4),'_ParMachine.mat']);
copyfile([pathname,filemot(1:end-4),'.mat'],[PathName,FileName(1:end-4),'.mat']);

