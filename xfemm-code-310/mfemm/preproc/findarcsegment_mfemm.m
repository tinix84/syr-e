function [ind] = findarcsegment_mfemm(FemmProblem, loc)
% finds the arc nearest to a given location loc=[x y];
%
% Syntax
%
% [ind] = findarcsegment_mfemm(FemmProblem, loc)
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

if ~isfield(FemmProblem,'ArcSegments')
    ind=[];
    return
end
if ~isfield(FemmProblem,'Nodes')
    ind=[];
    return
end
nodelist=FemmProblem.Nodes;
arclist=FemmProblem.ArcSegments;
if(numel(arclist)==0)
    ind=[];
    return
end
ind=1;
d0=shortestDistanceFromArc(FemmProblem,arclist(1),loc(1)+1i*loc(2));
for i=1:numel(arclist)
    d1=shortestDistanceFromArc(FemmProblem,arclist(i),loc(1)+1i*loc(2));
    if(d1<d0)
        d0=d1;
        ind=i;
    end
end