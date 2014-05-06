function [p,j]=getArcArcIntersection(FemmProblem, arc0, arc1)
% Returns the complex coordinates of the intersection(s) between two given
% arcs. Both arcs must be given in the form
% of FemmProblem fileds. The j value returned is the numeber of found
% intersections. It can be 0,1,2. Posiotion(s) of the intersection are
% returned as an array of complex numebers p(i)=x(i)+1i*y(i);
%
% Syntax
%
% [p,j]=getArcArcIntersection(FemmProblem, arc0, arc1)

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

j=0; p=[];
nodelist = FemmProblem.Nodes;

x=nodelist(arc0.n0+1).Coords;
a0=x(1)+1i*x(2);
x=nodelist(arc1.n0+1).Coords;
a1=x(1)+1i*x(2);

[c1,R1]=getCircle(FemmProblem,arc1);
[c0,R0]=getCircle(FemmProblem,arc0);

d=abs(c1-c0); %distance between centers
inds = d >= abs(R1 - R0) && d <= (R1 + R0);
if (sum(inds))==0
    return
end
if ((d>R0+R1) || (d<1.e-08))
    return
end
% directly eliminate case where there can't
% be any crossings....

l=sqrt((R0+R1-d)*(d+R0-R1)*(d-R0+R1)*(d+R0+R1))/(2.*d);
c=1.+(R0/d)*(R0/d)-(R1/d)*(R1/d);
t=(c1-c0)/d;
tta0=arc0.ArcLength*pi/180;
tta1=arc1.ArcLength*pi/180;

p(j+1)=c0 + (c*d/2.+ 1i*l)*t;		% first possible intersection;
z0=angle((p(j+1)-c0)/(a0-c0));
z1=angle((p(j+1)-c1)/(a1-c1));
if ((z0>0.) && (z0<tta0) && (z1>0.) && (z1<tta1)),j=j+1;end

if(abs(d-R0+R1)/(R0+R1)< 1.e-05)
    p(j+1)=c0+ c*d*t/2.;
    return
end

p(j+1)=c0 + (c*d/2.-1i*l)*t;		% second possible intersection
z0=angle((p(j+1)-c0)/(a0-c0));
z1=angle((p(j+1)-c1)/(a1-c1));
if ((z0>0.) && (z0<tta0) && (z1>0.) && (z1<tta1)),j=j+1;end

% returns the number of valid intersections found;
% intersections are returned in the array p[];
