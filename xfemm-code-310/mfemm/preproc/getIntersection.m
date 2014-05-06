function [xi,yi]=getIntersection(FemmProblem,n0,n1,segmentIndex)
% Returns the [xi,yi] coordinates of an intersection among the
% segmentIndex-th segment in FemmProblem and the segment fron node n0 to
% node n1.
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

linelist = FemmProblem.Segments;
nodelist = FemmProblem.Nodes;
%Check to see if the two lines have a common endpoint
%If they do, there can be no other intersection...
if (n0==linelist(segmentIndex).n0),
    xi=FemmProblem.Nodes(n0).Coords(1);
    yi=FemmProblem.Nodes(n0).Coords(2);
    return
end
if (n0==linelist(segmentIndex).n1)
    xi=FemmProblem.Nodes(n1).Coords(1);
    yi=FemmProblem.Nodes(n1).Coords(2);
    return
end
if (n1==linelist(segmentIndex).n0)
    xi=FemmProblem.Nodes(n0).Coords(1);
    yi=FemmProblem.Nodes(n0).Coords(2);
    return
end
if (n1==linelist(segmentIndex).n1)
    xi=FemmProblem.Nodes(n1).Coords(1);
    yi=FemmProblem.Nodes(n1).Coords(2);
    return
end

x=nodelist(linelist(segmentIndex).n0+1).Coords;
p0=x(1)+1i*x(2);
x=nodelist(linelist(segmentIndex).n1+1).Coords;
p1=x(1)+1i*x(2);
x=nodelist(n0+1).Coords;
q0=x(1)+1i*x(2);
x=nodelist(n1+1).Coords;
q1=x(1)+1i*x(2);
ee=min([abs(p1-p0),abs(q1-q0)])*1.0e-8;

% Rotate and scale the prospective line
q0=(q0-p0)/(p1-p0);
q1=(q1-p0)/(p1-p0);

% Check for cases where there is obviously no intersection
if ((real(q0)<=0) && (real(q1)<=0.))
    xi=[];yi=[];
    return
end
if ((real(q0)>=1.) && (real(q1)>=1.))
    xi=[];yi=[];
    return
end
if ((imag(q0)<=0.) && (imag(q1)<=0.))
    xi=[];yi=[];
    return
end
if ((imag(q0)>=0.) && (imag(q1)>=0.))
    xi=[];yi=[];
    return
end

% compute intersection
z=imag(q0)/imag(q0-q1);
%check to see if the line segments intersect at a point sufficiently
%far from the segment endpoints....
x=real((1.0 - z)*q0 + z*q1);
if((x < ee) || (x > (1.0 - ee)))
    xi=[];yi=[];
    return
end

% return resulting intersection point
p0 = (1.0 - z)*nodelist(n0+1).Coords + z*nodelist(n1+1).Coords;
xi=p0(1);
yi=p0(2);

