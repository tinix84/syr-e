function d=shortestDistanceFromArc(FemmProblem, arc, p)
% Returns the d minimum distance of a point p with complex coordinates 
% p=x+1i*y from arc. In this function arc is an arc field of FemmProblem.
% This is a Matlab porting of the equivalent cpp FEMM function.
%
% Syntax
%
% d=shortestDistanceFromArc(FemmProblem, arc, p)

% See also: 
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
x=nodelist(arc.n0+1).Coords;
a0=x(1)+1i*x(2);
x=nodelist(arc.n1+1).Coords;
a1=x(1)+1i*x(2);
[c,R]=getCircle(FemmProblem,arc);
d=abs(p-c);
if(d==0)
    d=R;
    return
end
t=(p-c)/d;
l=abs(p-c-R*t);
z=angle(t/(a0-c))*180/pi;
if((z>0)&&(z<arc.ArcLength))
    d=l;
    return
end
z=abs(p-a0);
l=abs(p-a1);
if(z<l)
    d=z;
    return
end
d=l;