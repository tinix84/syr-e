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

function dataSet = DrawPushMachine(handles,filename,pathname)

dataSet = handles.dataSet;

if nargin < 2
    [filename,pathname] = uiputfile(['newmachine.fem'],'input machine name and location');
else
    filename = strrep(filename,'.mat','.fem');
end

%% ==== FIRST PART FROM FEMMFitnessX ======================================
[~, ~, geo, per, mat] = data0(dataSet);
RQ = dataSet.RQ;
currentDir = pwd;

% RQ defines the candidate machine 
[geo,gamma,mat] = interpretRQ(RQ,geo,mat);

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
eval_type = 'singt';

% FEMM
openfemm
[geo,mat] = draw_motor_in_FEMM(geo,eval_type,mat);
mi_saveas([pathname filename]);
mi_close, closefemm

geo.RQ = RQ;

filename = strrep(filename,'fem','mat');
save([pathname filename],'FemmProblem','geo','per','dataSet','mat');

set(handles.currentMotFileName,'String',filename);  % update display

%% refresh GUI display data
load([pathname filename]);
dataSet.RQ = round(dataSet.RQ,4);
dataSet.currentpathname = [pathname '\'];
dataSet.currentfilename = filename;

cd(currentDir);

