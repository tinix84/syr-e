function [c,R]=getCircle(FemmProblem,arc)
% Returns the center -as a complex variable- and the radius of circle
% containg the given arc. In this function arc is an arc field of FemmProblem.
% This is a Matlab porting of the equivalent cpp FEMM function.
%
% Syntax
%
% [c,R]=getCircle(FemmProblem,arc)

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
d=abs(a1-a0);			% distance between arc endpoints

% figure out what the radius of the circle is...
t=(a1-a0)/d;
tta=arc.ArcLength*pi/180.;
R=d/(2.*sin(tta/2.));
c=a0 + (d/2. + 1i*sqrt(R*R-d*d/4.))*t; % center of the arc segment's circle...
