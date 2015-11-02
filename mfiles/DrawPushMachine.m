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

%% == NEW SCRIPT FOR SAVING A MACHINE FROM GUI MANUAL ENTRIES ==

%% ==== TAKE THE DIR AND THE STRUCT DATASET ===============================

dataSet = handles.dataSet;

[filename,dname] = uiputfile('newmachine.fem','input machine name and location');
dname = dname(1:end-1);

%% ==== FIRST PART FROM FEMMFitnessX ======================================

[bounds, geo, per] = data0(dataSet);
RQ = dataSet.RQ;
currentDir = pwd;

% RQ defines the candidate machine 
geo.pathname = cd;

[geo,gamma] = interpretRQ(RQ,geo);

FemmProblem = loadfemmfile([currentDir filesep 'empty_case.fem']);

FemmProblem.ProbInfo.Depth = geo.l;
FemmProblem.Segments = [];
FemmProblem.ArcSegments = [];
FemmProblem.Nodes = [];
FemmProblem.BoundaryProps = [];
FemmProblem.Circuits = [];
FemmProblem.BlockLabels = [];
FemmProblem.PointProps = [];
eval_type = 'MO_OA';
cd
% [geo,FemmProblem] = draw_motor_in_XFEMM(geo,eval_type,FemmProblem);
openfemm
[geo] = draw_motor_in_FEMM(geo,eval_type);
delete('mot0.fem');

%% SAVE THE FILE ==========================================================
cd(dname);
% writefemmfile(filename, FemmProblem);
mi_saveas([dname '\' filename]);
mi_close, closefemm
geo.RQ = RQ;
dataSet.RQnames = geo.RQnames;
dataSet.Dalpha1BouCheck = 0;
dataSet.DalphaBouCheck = 0;
dataSet.hcBouCheck = 0;
dataSet.DxBouCheck = 0;
dataSet.GammaBouCheck = 0;
dataSet.GapBouCheck  = 0;
dataSet.BrBouCheck  = 0;
dataSet.AirgapRadiusBouCheck  = 0;
dataSet.ToothWidthBouCheck  = 0;
dataSet.ToothLengthBouCheck  = 0;
dataSet.StatorSlotOpenBouCheck  = 0;
dataSet.ToothTangDepthBouCheck  = 0;
save(strrep(filename,'fem','mat'),'FemmProblem','geo','per','dataSet');
cd(currentDir);

