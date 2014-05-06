function [p,j]=getLineArcIntersection(FemmProblem, arc, segment)
% Returns the complex coordinates of the intersection(s) between a given
% segment and a given arc. Both arc and segment must be given in the form
% of FemmProblem fileds. The j value returned is the numeber of found
% intersections. It can be 0,1,2.
%
% Syntax
%
% [p,j]=getLineArcIntersection(FemmProblem, arc, segment)

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

j=0;
nodelist = FemmProblem.Nodes;
x=nodelist(segment.n0+1).Coords;
p0=x(1)+1i*x(2);
x=nodelist(segment.n1+1).Coords;
p1=x(1)+1i*x(2);
x=nodelist(arc.n0+1).Coords;
a0=x(1)+1i*x(2);
x=nodelist(arc.n1+1).Coords;
a1=x(1)+1i*x(2);
d=abs(a1-a0); %Distance between arc endpoints

% figure out what the radius of the circle is...
t=(a1-a0)/d;
tta=arc.ArcLength*pi/180;
R=d/(2*sin(tta/2));
c=a0 + (d/2. + 1i*sqrt(R*R-d*d/4.))*t; % center of the arc segment's circle...

% figure out the distance between line and circle's center;
d=abs(p1-p0);
t=(p1-p0)/d;
v=(c-p0)/t;
if (abs(imag(v))>R)
    j=0;
    p=[];
    return
end
l=sqrt( R*R - imag(v)*imag(v));	% imag(v) is distance between line and center...
if ((l/R) < 1.e-05)             % case where line is very close to a tangent;
    p(j+1)=p0 + real(v)*t;          % make it be a tangent.
    R=real((p(j+1)-p0)/t);
    z=angle((p(j+1)-c)/(a0-c));
    if ((R>0) && (R<d) && (z>0.) && (z<tta)) 
        j=j+1;
        return
    end
end

p(j+1)=p0 + (real(v)+l)*t;		% first possible intersection;
R=real((p(j+1)-p0)/t);
z=angle((p(j+1)-c)/(a0-c));
if ((R>0) && (R<d) && (z>0.) && (z<tta))
    j=j+1;
end

p(j+1)=p0 + (real(v)-l)*t;		% second possible intersection
R=real((p(j+1)-p0)/t);
z=angle((p(j+1)-c)/(a0-c));
if ((R>0) && (R<d) && (z>0.) && (z<tta)) 
    j=j+1;
end                              
% returns the number of valid intersections found;
% intersections are returned in the array p;
if j==0
    p=[];
end
     