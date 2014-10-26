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

%% 2013/08/21 MG script che prepara una macchina elaborata da MOGA per le successive simulazioni femm, Magnet ecc
% input: in questa fase ancora dati0 e il file.mat in risultati del motore
% di cui si vuole procedere alle successive rielaborazioni.
% Output: viene creata una cartella col nome del motore e 2 sottocartelle:
% DXF: contenete il file.mat e .dxf del motore...
% FEMM: contenente il file femm e .mat del motore...

clear all; close all; clc;

current_path=cd;
[pathstr, name, ext] = fileparts(current_path);

% addpath C:\Programmi\femm42\mfiles\;
% addpath('C:\Program Files (x86)\femm42\mfiles');
addpath('d:\femm42_beta\mfiles\');
addpath([pathstr,'\dxf_conv_fun\']);
addpath([pathstr,'\mfiles\']);

load ultimo.mat;

[filemot, pathname] = uigetfile([pathname '\*m*.mat'], 'Pick a motor');
save ultimo.mat pathname;
run([pathname,'data0']);
load([pathname filemot]);


eval_type='singt';
error_code = 0;
fem=dimMesh(geo,eval_type);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qs=geo.Qs;
ps=geo.ps;

if (size(geo.avv,2)<Qs)
    avv=[geo.avv,-geo.avv];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  %%%%%%
%% Viene salvato il file nel formato dxf
%%  %%%%%%
pathname_DXF=[pathname,filemot(1:end-4),'\DXF\'];
Export_to_DXF;

%%  %%%
%% Salvataggio del file.mat
%%  %%%
BLKLABELS;
Mac.g=geo.g;
Mac.R=geo.R;
Mac.RtS=geo.r+geo.g;
Mac.p=geo.p;
Mac.Q=geo.ns*geo.p;
Mac.q=Mac.Q/(6*geo.p);
Mac.l=geo.l;
Mac.Qs=Qs;
Mac.ps=ps;
Mac.th0=geo.th0;
Mac.N_turn=geo.Ns;
Mac.N_cond=geo.Ns/geo.p/Mac.q/size(avv,1); % numero di conduttori in cava per strato
Mac.N_parallel=1;
Mac.avv=avv;
Mac.n_mag_simulati=0;
Mac.fem=fem;
Mac.caso=1;
Mac.MachineName=filemot(1:end-4);
Mac.MachineNameMn=[filemot(1:end-4),'.mn'];
save([pathname_DXF,filemot(1:end-4),'_ParMachine.mat'],'Mac','BLKLABELS','statore','rotore2');
copyfile([pathname,filemot(1:end-4),'.mat'],[pathname_DXF,filemot(1:end-4),'.mat']);

%% %%%%%
%% Viene salvato il file nel formato femm
%% %%%%%
% pathname_FEMM=[pathname,filemot(1:end-4),'\FEMM\'];
% Export_to_FEMM;
% save([pathname_FEMM,filemot(1:end-4),'_ParMachine.mat'],'Mac','fem','BLKLABELS','statore','rotore','STATOREdati');

