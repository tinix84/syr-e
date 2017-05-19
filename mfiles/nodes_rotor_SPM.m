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

function [geo,mat,temp] = nodes_rotor_SPM(geo,mat) %,mat,fem)

% parameters

r = geo.r;                      % rotor radius
p = geo.p;                      % pole pairs
phi = geo.phi/p;                % angle range of permanent magnet
lm = geo.lm;                    % the thickness of permant magnet
hc = 0;
seg = geo.dx;                   % the number of segments of magnet

PMregular = geo.BarFillFac;
if PMregular > 1
    PMregular =1;                      % limit to per unit
end
%% COMPLETE DESIGN (ARCO OUTDOOR AIR GAP AND LAYER 1) AND ASSIGNMENTS (MATERIALS AND BOUNDARY CONDITIONS)
x = r * cos(90/p * pi/180);
y = r * sin(90/p * pi/180);
%% DISEGNO LO STRATO 1 DEL TRAFERRO (VA DA (th_m0) A (th_m0+180/p))

hybrid = 0;
if hybrid == 1
    %% hybrid shape
    PM_angle = phi/4;
    Fe_angle = phi/4;
    xPMco = r; yPMco = 0;
    [xPMo,yPMo] = rot_point(xPMco,yPMco,PM_angle*pi/180);
    [xFeo,yFeo] = rot_point(xPMo,yPMo,Fe_angle*pi/180);
    
    [x5, y5] = rot_point(xPMco,yPMco,90/p*pi/180);
    
    xPMci = r-lm; yPMci = 0;
    [xPMi,yPMi] = rot_point(xPMci,yPMci,PM_angle*pi/180);
    [xFei,yFei] = rot_point(xPMi,yPMi,Fe_angle*pi/180);
    [x4, y4] = rot_point(xPMci,yPMci,90/p*pi/180);
else
    %%  regular or rounded shape
    xPMco = r; yPMco = 0;   

    xPMregular = r-lm + PMregular*lm; yPMregular = 0;
    
    [xPMo,yPMo] = rot_point(xPMregular,yPMregular,phi/2*pi/180);    % PM edge point
    [x5, y5] = rot_point(xPMco,yPMco,90/p*pi/180);                  % Air zone point   
    [x6, y6] = rot_point(xPMco,yPMco,phi/2*pi/180);                 % Air zone point
    
    x6= x5;
    y6= y5;
    
    xPMci = r-lm; yPMci = 0;            
    [xPMi,yPMi] = rot_point(xPMci,yPMci,phi/2*pi/180);              % PM edge point on the steel
    
    [x4, y4] = rot_point(xPMci,yPMci,90/p*pi/180);                  % Rotor edge point   
        %% chao 2017.01.09 use an arc to assume sinusoidal
    xArccenter = (xPMco + xPMo - (yPMo^2/(xPMco-xPMo)))/2;          % find arc center location
    yArccenter = 0;
end
%% segmentation point build by sine wave
% if seg ~=1    
%     NoSeg = floor(seg);
%     for jj = 1:floor(NoSeg/2)
%         rPMso(jj) = r- lm*(1-PMregular)+ lm*(1-PMregular) * sin(pi/seg*jj);
%         xPMso(jj) = rPMso(jj) * cos(phi/2*pi/180-phi/2*pi/180*2/seg*jj);
%         yPMso(jj) = rPMso(jj) * sin(phi/2*pi/180-phi/2*pi/180*2/seg*jj);
%         [xPMsi(jj),yPMsi(jj)] = rot_point(xPMi,yPMi,-phi/2*pi/180*2/seg*jj);
%         temp.xPMso(jj) = xPMso(jj);
%         temp.yPMso(jj) = yPMso(jj);
%         temp.xPMsi(jj) = xPMsi(jj);
%         temp.yPMsi(jj) = yPMsi(jj);
%     end 
% end

%% segmentation point build by circular wave
if seg~=1
    % a quarter pole is considered
    NoSeg = 1:floor(seg/2);
    SegAngle = phi/seg;                         % angle span of each segment
    % coodinate displacement to Arccenter to calculate y on outside shape of PM
    [xPMso,yPMso] = rot_point(xPMo-xArccenter,yPMo,-atan(yPMo/(xPMo-xArccenter))*2/seg*NoSeg);
    % coodinate recovery
    xPMso = xPMso + xArccenter;
    % get positions on inside shape of PM
    [xPMsi,yPMsi] = rot_point(xPMi,yPMi,-(NoSeg*SegAngle)*pi/180);
    temp.xPMso = xPMso;
    temp.yPMso = yPMso;
    temp.xPMsi = xPMsi;
    temp.yPMsi = yPMsi;
end
%% calculate PM area
area2 = phi/2/360*pi*r^2 - 0.5*xArccenter*xPMo;
area3 = atan(yPMo/(xPMo-xArccenter))/2*(xPMco-xArccenter)^2 - area2;
PMarea = 2*(area3 + phi/2/360*pi*xPMregular^2 - phi/2/360*pi*xPMci^2);
geo.PMarea = PMarea * 2*p*1e-6;
geo.PMmass = geo.PMarea*geo.l*mat.LayerMag.kgm3*1e-3;

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
if hybrid == 0
    temp.xArccenter = xArccenter;
    temp.yArccenter = yArccenter;
    temp.x6 = x6;
    temp.y6 = y6;
else
    temp.xFeo = xFeo;
    temp.xFei = xFei;
    temp.yFeo = yFeo;
    temp.yFei = yFei;
    geo.PM_angle = PM_angle;
    geo.Fe_angle = Fe_angle;
    geo.hybrid = 1;
end
geo.hybrid = hybrid;
geo.hc = hc;

