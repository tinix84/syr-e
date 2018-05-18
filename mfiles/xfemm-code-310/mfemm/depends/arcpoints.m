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

function [x, y] = arcpoints(A, B, angle, maxdeg)
% get the points on an arc for plotting

    % convert the angles to radians
    angle = deg2rad(angle);
    maxdeg = deg2rad(maxdeg);

    % get centre and radius of circles
    [centre, r] = circcentre(A, B, angle);
    
    % get starting angle of arcs
    tempA = A - centre;
    
    [starttheta, rho] = cart2pol(tempA(:,1), tempA(:,2));
    
    % get arc points
    npnts = max(3, ceil(angle ./ maxdeg));
    
    pnts = linspace(starttheta, starttheta + angle, npnts);
    
    [x, y] = pol2cart(pnts, repmat(rho, size(pnts)));
    
    x = x + centre(:,1);
    y = y + centre(:,2);
    
    
end