% Copyright 2018
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

function dataSet = DrawPushMachine_MN(handles,filename,pathname)

dataSet = handles.dataSet;

if nargin < 2
    [filename,pathname] = uiputfile(['newmachine.mat'],'input machine name and location');
%     pathname = pathname(1:end-1);
end

[~, ~, geo, per, mat] = data0(dataSet);

eval_type = 'singt';

h = OpenMagnet(1);  % 1 = visible, 0 = invisible
[geo,mat] = draw_motor_in_MN(geo,eval_type,mat,pathname,filename,h);
[h,f] = SaveDocumentMagnet(h,[pathname,filename(1:end-4),'.mn']);
CloseMagnet(h)

save([pathname filename],'geo','per','dataSet','mat'); 

set(handles.currentMotFileName,'String',filename);  % update display

% refresh GUI display data
load([pathname filename]);
dataSet.RQ = round(dataSet.RQ,4);
dataSet.currentpathname = [pathname '\'];
dataSet.currentfilename = filename;

