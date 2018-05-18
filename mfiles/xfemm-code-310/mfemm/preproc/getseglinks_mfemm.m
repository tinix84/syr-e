function links = getseglinks_mfemm(FemmProblem)
% gets the links between nodes from an mfemm FemmProblem structure
%
% Syntax
%
% links = getlinks_mfemm(FemmProblem)
%
% Input
%
%   FemmProblem - structure containing a field 'Segments' which is an n
%     element structure array with fields n0 and n1 representing indices of
%     the nodes making up the segment.
%
% Output
%
%   links - (n x 2) array containing pairs of indices of nodes making up
%     segments
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


if isfield(FemmProblem, 'Segments')
    if ~isempty(FemmProblem.Segments)
        links = [cell2mat({FemmProblem.Segments(:).n0}'), ...
            cell2mat({FemmProblem.Segments(:).n1}')];
    else
        links = [];
    end
else
    links = [];
end

end