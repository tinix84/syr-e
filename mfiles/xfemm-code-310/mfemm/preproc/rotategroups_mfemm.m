function FemmProblem = rotategroups_mfemm(FemmProblem, groupnos, alfa)
% rotates all nodes, segments, arc segments and block labels which are
% members of the specified group numbers
%
% Syntax
%
% FemmProblem = rotategroups_mfemm(FemmProblem, groupnos, alfa)
%
%

% Copyright 2012 Richard Crozier
% Copyright 2014 Ernesto Mininno
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

excludenodes = [];
includenodes_id = [];

for i = 1:numel(groupnos)
    if isfield(FemmProblem, 'Segments')
        for ii = 1:numel(FemmProblem.Segments)
            if FemmProblem.Segments(ii).InGroup == groupnos(i)
                includenodes_id = [includenodes_id,...
                    FemmProblem.Segments(ii).n0,...
                    FemmProblem.Segments(ii).n1];
            end
        end
    end
    
    if isfield(FemmProblem, 'ArcSegments')
        for ii = 1:numel(FemmProblem.ArcSegments)
            if FemmProblem.ArcSegments(ii).InGroup == groupnos(i)
                includenodes_id = [includenodes_id,...
                    FemmProblem.ArcSegments(ii).n0,...
                    FemmProblem.ArcSegments(ii).n1];
            end
        end
    end
    
    if isfield(FemmProblem, 'BlockLabels')
        for ii = 1:numel(FemmProblem.BlockLabels)
            if FemmProblem.BlockLabels(ii).InGroup == groupnos(i)
                FemmProblem.BlockLabels(ii).Coords = (rotationMatrix*(FemmProblem.BlockLabels(ii).Coords'))';
            end
        end
    end
end

if isfield(FemmProblem, 'Nodes')
    for i = 1:numel(groupnos)
        for ii = 1:numel(FemmProblem.Nodes)
            if FemmProblem.Nodes(ii).InGroup == groupnos(i) && ...
                    ~any(excludenodes == ii)
                includenodes_id = [includenodes_id, ii-1];
            end
        end
    end
end

includenodes_id = unique(includenodes_id);
FemmProblem = rotatetenodes_mfemm(FemmProblem, alfa, includenodes_id);
end