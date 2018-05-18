function [id, xycoords] = findsegment_mfemm(FemmProblem, loc)
% finds the segment with mid point nearest a given location loc=[x y];
%
% Syntax
%
% [id, xycoords] = findsegment_mfemm(FemmProblem, loc)
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

    segcoords = getsegmidpointcoords_mfemm(FemmProblem);
    
    % find the nearest node to the location
    % ipdm returns a structure with fields named 'rowindex',
    % 'columnindex', and 'distance'.
    result = ipdm(loc, segcoords, 'Result', 'Structure', 'Subset', 'NearestNeighbor');
    
    % get the indices of the nodes and subtract 1 to make zero based
    id = result.columnindex - 1;
    
    % return the actual coordinates of the segment mid-points
    xycoords = cat(1, segcoords(id+1,:));

end