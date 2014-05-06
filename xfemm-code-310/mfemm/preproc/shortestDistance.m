function d=shortestDistance(FemmProblem, p , q, segmentIndex)
% Returns the d minimum distance of a point with  coordinates 
% (p,q) from the segmentIndex-th segment in FemmProblem.
% This is a Matlab porting of the equivalent cpp FEMM function.
%
% Syntax
%
% [xi,yi]=getIntersection(FemmProblem,n0,n1,segmentIndex)

% See also: newarcsegment_mfemm.m
%
% Copyright 2012 Richard Crozier
% Copyright 2014 Ernesto Mininno (e.mininno@gmail.com)

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
nodelist=FemmProblem.Nodes;
linelist=FemmProblem.Segments;
x(1)=nodelist(linelist(segmentIndex).n0+1).Coords(1);
y(1)=nodelist(linelist(segmentIndex).n0+1).Coords(2);
x(2)=nodelist(linelist(segmentIndex).n1+1).Coords(1);
y(2)=nodelist(linelist(segmentIndex).n1+1).Coords(2);
t=((p-x(1))*(x(2)-x(1)) + (q-y(1))*(y(2)-y(1)))/...
    ((x(2)-x(1))*(x(2)-x(1)) + (y(2)-y(1))*(y(2)-y(1)));

if (t>1)
    t=1;
end
if (t<0)
    t=0;
end

x(3)=x(1)+t*(x(2)-x(1));
y(3)=y(1)+t*(y(2)-y(1));
d = sqrt((p-x(3))*(p-x(3)) + (q-y(3))*(q-y(3)));
end