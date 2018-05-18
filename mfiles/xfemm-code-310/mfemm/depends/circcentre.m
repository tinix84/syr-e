% Copyright 2014
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

function [centre, r] = circcentre(A, B, angle)
% calculates the centre and radius of a circle given two points and an arc
% angle between them. The position of the circle is determined by the
% order of the supplied points.
%
% Syntax
%
% [centre, r] = circcentre(A, B, angle)
%
%
    
    % get vector pointing from A to B
    AB = B - A;
    
    % find perpendicular vector to AB
    V = [ -AB(:,2), AB(:,1) ];
    
    % find mid point of AB
    M = A + AB .* 0.5;
    
    % find length of AB and divide by two to get triangle base
    b = 0.5 * magn(AB);
    
    % find triangle height
    h = b ./ tan(angle ./ 2);
    
    % find circle centre
    centre = M + h * unit(V);
    
    % find radius
    r = sqrt(h.^2 + b.^2);
    
end