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

%% 2013/07/05 MG costruzione della macchina in FEMM note le matrici di
%% statore e rotore
%% Si disegna e assegnano le proprietà secondo la struttura matriciale:

%% Definizione dei gruppi:
groupStat=1;
groupRot=2;
groupTraf=20;
closefemm;
filename = 'mot0.fem';
% opendocument('empty_case.fem',h_temp);
openfemm;
opendocument([pathstr,'\empty_case.fem']);

mi_probdef(0,'millimeters','planar',1e-8,Mac.l,25);

%% Scelta periodicità per i boundary:
% if (gcd(Mac.Qs,Mac.ps)>1)
%     Period=4;
% else
%     Period=5;
% end
Period=5;

% boundary
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


%Disegna rotore;
Importa_Matrice(rotore,groupRot,fem.res);
%Assegna materiali rotore;

Assegna_material_rotor(BLKLABELS,fem,groupRot);


% Disegna statore
Importa_Matrice(statore,groupStat,fem.res);
% Assegna materiali statore
Assegna_material_stat_1Nbob(BLKLABELS,Mac,fem,groupStat);

% Disegno traferro:
pc=360/(Mac.Q)/2; % mezzo passo cava di statore
xr=Mac.RtS-Mac.g;
AirGapBuild(Mac.Qs,Mac.ps,Mac.p,Mac.g,pc,xr,fem.res_traf,groupStat,groupRot,groupTraf);
disegna_bordo_mobile_polPlus(Mac,th_m0,fem,Mac.ps);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assegnazione condizioni al contorno
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% per il momento non è buono perchè
% sarebbe meglio una selezione automatica e invece?!?!....:-(
% Condizioni al conrtorno di statore

BLKLABELSstat=BLKLABELS.statore;
for ii=1:size(BLKLABELSstat.boundary,1)
    
    mi_selectsegment(BLKLABELSstat.boundary(ii,1),BLKLABELSstat.boundary(ii,2));
    if (BLKLABELSstat.boundary(ii,3)==10)
        mi_setsegmentprop('APs1', 0, 1, 0, groupStat);
        mi_clearselected;
    elseif(BLKLABELSstat.boundary(ii,3)==0)
        mi_selectarcsegment(BLKLABELSstat.boundary(ii,1),BLKLABELSstat.boundary(ii,2))
        mi_setarcsegmentprop(fem.res, 'A=0', 0, groupStat);
        mi_clearselected;
    end
end
% Condizioni al conrtorno di rotore albero
BLKLABELSrot=BLKLABELS.rotore;
for ii=1:2
    
    mi_selectsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2));
    if (BLKLABELSrot.boundary(ii,3)==10)
        mi_setsegmentprop('APr1', 0, 1, 0, groupRot);
        mi_clearselected;
    elseif(BLKLABELSrot.boundary(ii,3)==0)
        mi_selectarcsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2))
        mi_setarcsegmentprop(fem.res, 'A=0', 0, groupRot);
        mi_clearselected;
    end
end
% Condizioni al conrtorno di rotore ferro lamierino
for ii=3:4
    
    mi_selectsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2));
    if (BLKLABELSrot.boundary(ii,3)==10)
        mi_setsegmentprop('APr2', 0, 1, 0, groupRot);
        mi_clearselected;
    elseif(BLKLABELSrot.boundary(ii,3)==0)
        mi_selectarcsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2))
        mi_setarcsegmentprop(fem.res, 'A=0', 0, groupRot);
        mi_clearselected;
    end
end

mkdir(pathname_FEMM);
mi_saveas([pathname_FEMM,filemot(1:end-4),'.fem']);

mi_close;closefemm;
