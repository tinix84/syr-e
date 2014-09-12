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

function FemmProblem=assign_block_prop_rotX(FemmProblem,BLKLABELS,fem,group)
BLKLABELSrot=BLKLABELS.rotore;

% Assegna aria alle barriere di flux:
for kk=1:size(BLKLABELSrot.names.BarName,1)
    %mi_addblocklabel(BLKLABELSrot.xy(kk,1),BLKLABELSrot.xy(kk,2));
    %mi_selectlabel(BLKLABELSrot.xy(kk,1),BLKLABELSrot.xy(kk,2));
    %mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(kk,3)}, 0, fem.res,'None', 0, group, 0);
    %mi_clearselected;
    FemmProblem = addblocklabel_mfemm(FemmProblem, BLKLABELSrot.xy(kk,1),BLKLABELSrot.xy(kk,2), ...
                                  'BlockType', BLKLABELS.materials{BLKLABELSrot.xy(kk,3)}, ...
                                  'MaxArea', fem.res,...
                                  'InGroup', group,...
                                  'Turns',1);
end

index=kk+1;
% Assegna ferro di rotore
%     mi_addblocklabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
%     mi_selectlabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
%     mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(index,3)}, 0, fem.res,'None', 0, group, 0);
%     mi_clearselected;
    FemmProblem = addblocklabel_mfemm(FemmProblem, BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2), ...
                                  'BlockType', BLKLABELS.materials{BLKLABELSrot.xy(index,3)}, ...
                                  'MaxArea', fem.res,...
                                  'InGroup', group,...
                                  'Turns',0);
    
index=index+1;
% Assegna materiale albero
%     mi_addblocklabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
%     mi_selectlabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
%     mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(index,3)}, 0, fem.res,'None', 0, group, 0);
%     mi_clearselected;
    FemmProblem = addblocklabel_mfemm(FemmProblem, BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2), ...
                                  'BlockType', BLKLABELS.materials{BLKLABELSrot.xy(index,3)}, ...
                                  'MaxArea', fem.res,...
                                  'InGroup', group,...
                                  'Turns',0);

end