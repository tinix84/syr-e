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

function [rotor,BLKLABELSrot,geo,mat] = ROTmatr(geo,fem,mat)

% Rotor construction.
% rotor:                  	one row per FEMM line or arc
% BLKLABELSrot.xy:          center points of FEMM blocks
% BLKLABELSrot.boundary:    one row per FEMM bounday condition
% BLKLABELSrot.BarName:     names of flux barrier blocks

Ar = geo.Ar;
r = geo.r;
ps = geo.ps;
p = geo.p;
lm = geo.lm;

if ~strcmp(geo.RotType,'SPM')
    mat.LayerMag.Br = mat.LayerMag.Br.*ones(1,geo.nlay);   % replicate Br in case it is scalar
end

switch geo.RotType
    case 'Circular'
        % build nodes, lines and arcs for half a pole
        if (0)  % select 0 if you want to use the new circular geometry, 1 if you want use the old circular. Old circular not updated
            [geo,mat,temp] = nodes_rotor_Circ(geo,mat);
            rotor = build_matrix_Circ(temp,geo);
            disp('Old circular is selected')
        else
            [geo,mat,temp] = nodes_rotor_Circ_dx(geo,mat);
            rotor = build_matrix_Circ_dx(temp,geo);
        end
    case 'ISeg'
        % build nodes, lines and arcs for half a pole
        [geo,mat,temp] = nodes_rotor_ISeg(geo,mat);
        rotor = build_matrix_ISeg(temp,geo);
    case 'Seg'
        % build nodes, lines and arcs for half a pole
        [geo,mat,temp] = nodes_rotor_Seg(geo,mat);
        rotor = build_matrix_Seg(temp,geo);
    case 'Fluid'
        % build nodes, lines and arcs for half a pole
        [geo,mat,temp] = nodes_rotor_Fluid(geo,mat);
        rotor = build_matrix_Fluid(temp,geo);
    case 'SPM'
        % build nodes, lines and arcs for half a pole
        [geo,mat,temp] = nodes_rotor_SPM(geo,mat);
        rotor = build_matrix_SPM(temp,geo);
end
% replicate the half pole
rotor = finishRotorMatrix(rotor,geo);

% find the centers of all blocks
BarCenter = defineBlockCenters(temp,fem,geo);
% Assign label names
BarName = defineBlockNames(temp,geo,mat);
% BOUNDARY CONDITIONS
if (ps<2*geo.p)
    codBound_periodic = 10;           % 10 = Odd or Even Periodicity
else
    codBound_periodic = -10;          % -10 = no periodicity, simulate full machine
end

% shaft boundary
[xShaftBound1,yShaftBound1] = rot_point(mean([0,Ar]),0,-90/p*pi/180);
[xShaftBound2,yShaftBound2] = rot_point(mean([0,Ar]),0,(ps-1/2)*180/p*pi/180);
% rotor boundary
if strcmp(geo.RotType,'SPM')
    [xRotBound1,yRotBound1] = rot_point(mean([Ar,r]),0,-90/p*pi/180);
    [xRotBound2,yRotBound2] = rot_point(mean([Ar,r]),0,(ps-1/2)*180/p*pi/180);
else
    [xRotBound1,yRotBound1] = rot_point(mean([Ar,r-lm]),0,-90/p*pi/180);          % for SPM motor
    [xRotBound2,yRotBound2] = rot_point(mean([Ar,r-lm]),0,(ps-1/2)*180/p*pi/180);
end


%%% OUTPUT DATA %%%
%%%%%%%%%%%%%%%%%%%

%%% rotore2 %%%
% Rotate rotor in zero position
[nrig,ncol] = size(rotor);
rotor2=[];
for ii=1:2:ncol-2
    [xtemp,ytemp]=rot_point(rotor(:,ii),rotor(:,ii+1),90/p*pi/180);
    rotor2=[rotor2,xtemp,ytemp];
end
rotor2=[rotor2,rotor(:,ncol)];
rotor = rotor2;
rotor=checkPlotMatrix(rotor,1e-9);

%%% Block centers %%%
BLKLABELSrot.xy     =   BarCenter;
BLKLABELSrot.BarName =   BarName';
% Rotate block labels selection points
[xtemp,ytemp]=rot_point(BLKLABELSrot.xy(:,1),BLKLABELSrot.xy(:,2),90/p*pi/180);
BLKLABELSrot.xy=[xtemp,ytemp,BLKLABELSrot.xy(:,3:end)];
clear xtemp ytemp;

% Magnetization direction rotation
xtemp=cos(atan2(BLKLABELSrot.xy(:,7),BLKLABELSrot.xy(:,6))+(pi/2/p-eps));
ytemp=sin(atan2(BLKLABELSrot.xy(:,7),BLKLABELSrot.xy(:,6))+(pi/2/p-eps));
BLKLABELSrot.xy(:,6)=xtemp;
BLKLABELSrot.xy(:,7)=ytemp;
clear xtemp ytemp;

%%% Boundaries %%%
BLKLABELSrot.boundary = [xShaftBound1,yShaftBound1,codBound_periodic;
    xShaftBound2,yShaftBound2,codBound_periodic;
    xRotBound1,yRotBound1,codBound_periodic;
    xRotBound2,yRotBound2,codBound_periodic];
% Rotate boundary selection points
[xtemp,ytemp]=rot_point(BLKLABELSrot.boundary(:,1),BLKLABELSrot.boundary(:,2),90/p*pi/180);
BLKLABELSrot.boundary=[xtemp,ytemp,BLKLABELSrot.boundary(:,3:end)];
%%

