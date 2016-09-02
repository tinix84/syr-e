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

function dataSet = DrawPushMachine(handles,filename,dname)

dataSet = handles.dataSet;

if nargin < 2
    [filename,dname] = uiputfile(['newmachine.fem'],'input machine name and location');
    dname = dname(1:end-1);
else
    filename = strrep(filename,'.mat','.fem');
end

%% ==== FIRST PART FROM FEMMFitnessX ======================================
[~, ~, geo, per, mat] = data0(dataSet);
RQ = dataSet.RQ;
currentDir = pwd;

% RQ defines the candidate machine 
% geo.pathname = cd;

[geo,gamma,mat] = interpretRQ(RQ,geo,mat);

% FemmProblem = loadfemmfile([currentDir filesep 'empty_case.fem']);
FemmProblem.ProbInfo.Frequency = 0;
FemmProblem.ProbInfo.Precision = 1e-8;
FemmProblem.ProbInfo.MinAngle = 15;
FemmProblem.ProbInfo.LengthUnits = 'millimeters';
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
openfemm
[geo,mat] = draw_motor_in_FEMM(geo,eval_type,mat);
delete('mot0.fem');

%% SAVE THE FILE ==========================================================
cd(dname);
% writefemmfile(filename, FemmProblem);
mi_saveas([dname '\' filename]);
mi_close, closefemm
geo.RQ = RQ;
% dataSet.RQnames = geo.RQnames;
% dataSet.Dalpha1BouCheck = 0;
% dataSet.DalphaBouCheck = 0;
% dataSet.hcBouCheck = 0;
% dataSet.DxBouCheck = 0;
% dataSet.GammaBouCheck = 0;
% dataSet.GapBouCheck  = 0;
% dataSet.BrBouCheck  = 0;
% dataSet.AirgapRadiusBouCheck  = 0;
% dataSet.ToothWidthBouCheck  = 0;
% dataSet.ToothLengthBouCheck  = 0;
% dataSet.StatorSlotOpenBouCheck  = 0;
% dataSet.ToothTangDepthBouCheck  = 0;

filename = strrep(filename,'fem','mat');
save(filename,'FemmProblem','geo','per','dataSet','mat');

set(handles.currentMotFileName,'String',filename);  % update display

%% refresh GUI display data
load([dname '\' filename]);
dataSet.RQ = roundn(dataSet.RQ,-4);
dataSet.currentpathname = [dname '\'];
dataSet.currentfilename = filename;

cd(currentDir);

