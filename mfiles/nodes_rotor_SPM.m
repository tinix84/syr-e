% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function [geo,temp] = nodes_rotor_SPM(geo) %,mat,fem)

% group0 = 0;                      % define the rotor as group 0 in FEMM

% parameters

r = geo.r;                      % rotor radius
p = geo.p;                      % pole pairs
phi = geo.phi/p;                % angle range of permanent magnet
lm = geo.lm;                    % the thickness of permant magnet
hc = 0;
%% COMPLETE DESIGN (ARCO OUTDOOR AIR GAP AND LAYER 1) AND ASSIGNMENTS (MATERIALS AND BOUNDARY CONDITIONS)
x = r * cos(90/p * pi/180);
y = r * sin(90/p * pi/180);
%% DISEGNO LO STRATO 1 DEL TRAFERRO (VA DA (th_m0) A (th_m0+180/p))

xPMco = r; yPMco = 0;
[xPMo,yPMo] = rot_point(xPMco,yPMco,phi/2*pi/180);

[x5, y5] = rot_point(xPMco,yPMco,90/p*pi/180);

xPMci = r-lm; yPMci = 0;
[xPMi,yPMi] = rot_point(xPMci,yPMci,phi/2*pi/180);

[x4, y4] = rot_point(xPMci,yPMci,90/p*pi/180);

%%  SAVE THE FINAL DATA:

temp.xPMco = xPMco;
temp.yPMco = yPMco;
temp.xPMo = xPMo;
temp.yPMo = yPMo;
temp.xPMci = xPMci;
temp.yPMci = yPMci;
temp.xPMi = xPMi;
temp.yPMi = yPMi;
temp.x4 = x4;
temp.y4 = y4;
temp.x5 = x5;
temp.y5 = y5;

geo.hc = hc;

