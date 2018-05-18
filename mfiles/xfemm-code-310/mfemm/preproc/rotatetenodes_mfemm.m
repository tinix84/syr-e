function [FemmProblem] = rotatetenodes_mfemm(FemmProblem, alfa, nodeids)
% rotates all nodes which are members of the specified nodeids id list.
%
% Syntax
%
% [FemmProblem] = rotatetenodes_mfemm(FemmProblem, alfa, nodeids)
%
% 

% Copyright 2012 Richard Crozier
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


    rotationMatrix=[cos(alfa) -sin(alfa);sin(alfa) cos(alfa)];
    if nargin < 3
        nodeinds = 1:numel(FemmProblem.Nodes);
    else
        nodeinds = nodeids + 1;
    end
    
    for ind = nodeinds
        FemmProblem.Nodes(ind).Coords = (rotationMatrix*FemmProblem.Nodes(ind).Coords')';
    end

end