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
geo.delta_FBS=0; % no pole deformation
flagCirc=0; % if 1, use the old circular geometry, without dx
fem = dimMesh(geo,'singt');
matFBS=mat;
geoFBS=geo;
if ~strcmp(geo.RotType,'SPM')
    mat.LayerMag.Br = mat.LayerMag.Br.*ones(1,geo.nlay);   % replicate Br in case it is scalar
end
%% 1) Find the design points (nodes_rotor_xxx)
switch geo.RotType
    case 'Circular'
        % build nodes, lines and arcs for half a pole
        if (flagCirc)  % select 0 if you want to use the new circular geometry, 1 if you want use the old circular. Old circular not updated
            [geo,mat,temp] = nodes_rotor_Circ(geo,mat);
            disp('Old circular is selected')
        else
            [geo,mat,temp] = nodes_rotor_Circ_dx(geo,mat);
        end
    case 'ISeg'
        % build nodes, lines and arcs for half a pole
        [geo,mat,temp] = nodes_rotor_ISeg(geo,mat);
    case 'Seg'
        % build nodes, lines and arcs for half a pole
        [geo,mat,temp] = nodes_rotor_Seg(geo,mat);
    case 'Fluid'
        % build nodes, lines and arcs for half a pole
        [geo,mat,temp] = nodes_rotor_Fluid(geo,mat);
    case 'SPM'
        % build nodes, lines and arcs for half a pole
        [geo,mat,temp] = nodes_rotor_SPM(geo,mat);
    case 'Vtype'
        % build nodes, lines and arcs for half a pole
        [geo,mat,temp] = nodes_rotor_Vtype(geo,mat);
end

%% 2) if no FBS, build rotor matrix, find BLKLABELS and BarNames
if geo.th_FBS==0
    switch geo.RotType
        case 'Circular'
            if flagCirc
                rotor=build_matrix_Circ(temp,geo);
            else
                rotor=build_matrix_Circ_dx(temp,geo);
            end
        case 'ISeg'
            rotor=build_matrix_ISeg(temp,geo);
        case 'Seg'
            rotor=build_matrix_Seg(temp,geo);
        case 'Fluid'
            rotor=build_matrix_Fluid(temp,geo);
        case 'SPM'
            rotor=build_matrix_SPM(temp,geo);
        case 'Vtype'
            rotor=build_matrix_Vtype(temp,geo);
    end
    % replicate the half pole
    rotor = finishRotorMatrix(rotor,geo);
    % find the centers of all blocks
    BarCenter = defineBlockCenters(temp,fem,geo);
else
    delta_FBS=geo.th_FBS*[-1 1];
    %matFBS=mat;
    rotor=[];
    BarCenter=[];
    for ii=1:length(delta_FBS)
        % Design of one deformed pole
        geoFBS=geo;
        geoFBS.ps=1;
        geoFBS.delta_FBS=delta_FBS(ii);
        [rotorFBS,BLKLABELSrotFBS,]=drawPoleFBS(geoFBS,matFBS,delta_FBS(ii));
        % pole rotation angle
        th_rot=pi/geo.p*(ii-1)+sum(delta_FBS(1:ii-1))+delta_FBS(ii)/2;
        % geometry rotation
        rotorTmp=[];
        for kk=1:2:size(rotorFBS,2)-2
            [xtemp,ytemp]=rot_point(rotorFBS(:,kk),rotorFBS(:,kk+1),th_rot);
            rotorTmp=[rotorTmp,xtemp,ytemp];
        end
        rotorTmp=[rotorTmp,rotorFBS(:,end)];
        rotorFBS=rotorTmp;
        % labels rotation
        xyLabels=BLKLABELSrotFBS.xy;
        for kk=1:length(xyLabels(:,1))
            [xtemp,ytemp]=rot_point(xyLabels(kk,1),xyLabels(kk,2),th_rot);
            th_mag=atan2(xyLabels(kk,7),xyLabels(kk,6));
            th_mag=th_mag+th_rot+pi*floor(rem(ii+1,2));
            xmag=cos(th_mag);
            ymag=sin(th_mag);
            xyLabels(kk,:)=[xtemp,ytemp,xyLabels(kk,3),xyLabels(kk,4),xyLabels(kk,5),xmag,ymag,xyLabels(kk,end)];
        end
        % add pole to the total geometry
        rotor=[rotor;rotorFBS];
        if ii==length(delta_FBS) % the last two rows of xyLabels contain the rotor iron and shaft points
            BarCenter=[BarCenter;xyLabels];
        else
            BarCenter=[BarCenter;xyLabels(1:end-2,:)];
        end
        % Check the rotor matrix, to avoid plot error due to arcs
        [rotor]=checkPlotMatrix(rotor,1e-9);
        % Add shaft and rotor outer radius
        rotor=[rotor;
            0,0,geo.r,0,geo.r*cos(2*pi/geo.p),geo.r*sin(2*pi/geo.p),1;
            0,0,geo.Ar,0,geo.Ar*cos(2*pi/geo.p),geo.Ar*sin(2*pi/geo.p),1;
            0,0,geo.r,0,NaN,NaN,0;
            0,0,geo.r*cos(2*pi/geo.p),geo.r*sin(2*pi/geo.p),NaN,NaN,0];
    end
end
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
if geo.th_FBS==0
    [nrig,ncol] = size(rotor);
    rotor2=[];
    for ii=1:2:ncol-2
        [xtemp,ytemp]=rot_point(rotor(:,ii),rotor(:,ii+1),90/p*pi/180);
        rotor2=[rotor2,xtemp,ytemp];
    end
    rotor2=[rotor2,rotor(:,ncol)];
    rotor = rotor2;
    [xtemp,ytemp]=rot_point(BarCenter(:,1),BarCenter(:,2),90/p*pi/180);
    BarCenter=[xtemp,ytemp,BarCenter(:,3:end)];
    % Magnetization direction rotation
    xtemp=cos(atan2(BarCenter(:,7),BarCenter(:,6))+(pi/2/p-eps));
    ytemp=sin(atan2(BarCenter(:,7),BarCenter(:,6))+(pi/2/p-eps));
    BarCenter(:,6)=xtemp;
    BarCenter(:,7)=ytemp;
    clear xtemp ytemp;
end
%%% Block centers %%%
BLKLABELSrot.xy     =   BarCenter;
BLKLABELSrot.BarName =   BarName';
% Rotate block labels selection points
% [xtemp,ytemp]=rot_point(BLKLABELSrot.xy(:,1),BLKLABELSrot.xy(:,2),90/p*pi/180);
% BLKLABELSrot.xy=[xtemp,ytemp,BLKLABELSrot.xy(:,3:end)];
% clear xtemp ytemp;



%%% Boundaries %%%
BLKLABELSrot.boundary = [xShaftBound1,yShaftBound1,codBound_periodic;
    xShaftBound2,yShaftBound2,codBound_periodic;
    xRotBound1,yRotBound1,codBound_periodic;
    xRotBound2,yRotBound2,codBound_periodic];
% Rotate boundary selection points
[xtemp,ytemp]=rot_point(BLKLABELSrot.boundary(:,1),BLKLABELSrot.boundary(:,2),90/p*pi/180);
BLKLABELSrot.boundary=[xtemp,ytemp,BLKLABELSrot.boundary(:,3:end)];
%%

