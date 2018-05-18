function [FemmProblem, seginds, nodeinds, nodeids] = addpolygon_mfemm(FemmProblem, coords, varargin)
% adds a polygon to an mfemm FemmProblem structure   
%
% Syntax
%
% [FemmProblem, seginds, nodeinds, nodeids] = addpolygon_mfemm(FemmProblem, coords, 'Parameter', 'Value')
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

    if size(coords, 2) > 2
        
        error('MFEMM:addpolygon_mfemm:badcoords', ...
            'Coordinates should be an n x 2 matrix.')
        
    elseif size(coords, 1) < 3
        
        error('MFEMM:addpolygon_mfemm:noenoughnodes', ...
            'You must specify at least three nodes to make a polygon.')
    else
        
        % add the polygon!
        
        % first add the necessary nodes at the vertices
        [FemmProblem, nodeinds, nodeids] = addnodes_mfemm(FemmProblem, coords(:,1), coords(:,2));
        
        % next create segments linking each new node to the following new
        % node, and then linking the last new node with the first new node
        [FemmProblem, seginds] = addsegments_mfemm(FemmProblem, ...
                                                   nodeids, ...
                                                   [nodeids(2:end), nodeids(1)], ...
                                                   varargin{:});
        
    end


end