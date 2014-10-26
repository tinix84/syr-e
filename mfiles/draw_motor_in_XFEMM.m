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

function [geo,FemmProblem] = draw_motor_in_XFEMM(geo,eval_type,FemmProblem)

% draw_motor_in_XFEMM.m
% builds the .fem motor model with the rotor in position zero

% input: see function argument
% output:
% - updated geo (individual for the machine)
% - mot0.fem (th_m = 0, i123 = [0,0,0])

fem = dimMesh(geo,eval_type);
fem.res_traf=0.5;
fem.res=4;

% machine periodicity (t) and number of poles to be simulated (ps)
Q=geo.ns*geo.p;
t=gcd(round(geo.ns*geo.p),geo.p);  % periodicity
if ((6*t/Q)>1)
    ps=2*geo.p/t;   % periodic machine
    Qs=Q/t;
    tempWinTable = geo.avv;
    periodicity=4;
else
    ps=geo.p/t;     % anti-periodic machine
    Qs=Q/2/t;
    tempWinTable = [geo.avv -geo.avv];
    periodicity=5;
end
geo.Qs=Qs;  % # of simulated slots
geo.ps=ps;  % # of simulated poles

% calc winding factor (kavv) and rotor offset (phase1_offset)
[kavv, phase1_offset] = calcKwTh0(tempWinTable);

% offset angle for coordinate transformations
th_m0 = 0;                              % rotor position [mec deg]
geo.th0 = th_m0*geo.p - phase1_offset;  % d- to alpha-axis offset [elt deg]

% Boundary conditions
% if (rem(geo.ps,2)==0)
%     periodicity=4;
% else
%     periodicity=5;
% end

%% definition of the boundary conditions
% inner and outer circles
[FemmProblem, boundindA0, boundnameA0]=addboundaryprop_mfemm(FemmProblem,'A=0',0);
% Periodicity or Anti-Periodicity (2 x rotor + 2 x stator + 3 x airgap + 1 x APmove) APmove is the sliding contour
[FemmProblem, boundindAPr1, boundnameAPr1]=addboundaryprop_mfemm(FemmProblem,'APr1',periodicity);
[FemmProblem, boundindAPr2, boundnameAPr2]=addboundaryprop_mfemm(FemmProblem,'APr2',periodicity);
[FemmProblem, boundindAPg1, boundnameAPg1]=addboundaryprop_mfemm(FemmProblem,'APg1',periodicity);
[FemmProblem, boundindAPg2, boundnameAPg2]=addboundaryprop_mfemm(FemmProblem,'APg2',periodicity);
[FemmProblem, boundindAPg3, boundnameAPg3]=addboundaryprop_mfemm(FemmProblem,'APg3',periodicity);
[FemmProblem, boundindAPmove, boundnameAPmove]=addboundaryprop_mfemm(FemmProblem,'APmove',periodicity);
[FemmProblem, boundindAPs1, boundnameAPs1]=addboundaryprop_mfemm(FemmProblem,'APs1',periodicity);
[FemmProblem, boundindAPs2, boundnameAPs2]=addboundaryprop_mfemm(FemmProblem,'APs2',periodicity);

%% ROTOR
% build the matrices which describe the rotor
geo.x0 = geo.r/cos(pi/2/geo.p);
ROTmatr;
BLKLABELS.rotore=BLKLABELSrot;
% draw lines and arcs
FemmProblem=draw_lines_archesX(FemmProblem,rotore2,2,fem.res);
% assign block labels
FemmProblem=assign_block_prop_rotX(FemmProblem,BLKLABELS,fem,2);
% assign boundary conditions (??)

for ii=1:2
    if (BLKLABELSrot.boundary(ii,3)==10)
        [id, xycoords] = findsegment_mfemm(FemmProblem, [BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2)]);
        FemmProblem.Segments(id+1).BoundaryMarker = boundnameAPr1;
    elseif(BLKLABELSrot.boundary(ii,3)==0)
        [id] = findarcsegment_mfemm(FemmProblem, [BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2)]);
        FemmProblem.ArcSegments(id+1).BoundaryMarker = boundnameA0;
    end
end
for ii=3:4
    if (BLKLABELSrot.boundary(ii,3)==10)
        [id, xycoords] = findsegment_mfemm(FemmProblem, [BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2)]);
        FemmProblem.Segments(id+1).BoundaryMarker = boundnameAPr2;
    elseif(BLKLABELSrot.boundary(ii,3)==0)
        [id] = findarcsegment_mfemm(FemmProblem, [BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2)]);
        FemmProblem.ArcSegments(id+1).BoundaryMarker = boundnameA0;
    end
end

%% STATOR
% builds the matrixes which describe the stator
STATmatr;
BLKLABELS.statore=BLKLABELSstat;
% draw lines and arcs
FemmProblem=draw_lines_archesX(FemmProblem,statore,1,fem.res);
% assign block labels
FemmProblem=assign_block_prop_statX(FemmProblem,BLKLABELS,geo,fem,1); % assegna materiali;
% assign boundary conditions
BLKLABELSstat=BLKLABELS.statore;
for ii=1:size(BLKLABELSstat.boundary,1)
    [id, xycoords] = findsegment_mfemm(FemmProblem, [BLKLABELSstat.boundary(ii,1),BLKLABELSstat.boundary(ii,2)]);
    if (BLKLABELSstat.boundary(ii,3)==10)
        FemmProblem.Segments(id+1).BoundaryMarker = boundnameAPs1;
    elseif(BLKLABELSstat.boundary(ii,3)==0)
        [ind] = findarcsegment_mfemm(FemmProblem, [BLKLABELSstat.boundary(ii,1),BLKLABELSstat.boundary(ii,2)]);
        FemmProblem.ArcSegments(ind).BoundaryMarker = boundnameA0;
    end
end

%% airgap (group 20)
    FemmProblem=AirGapBuildX(FemmProblem,Qs,ps,geo.p,geo.g,360/(geo.ns*geo.p)/2,geo.r,fem.res_traf,1,2,boundnameAPg1,boundnameAPg2,boundnameAPg3);
    FemmProblem=draw_airgap_arc_with_meshX(FemmProblem,geo,th_m0,fem,boundnameAPmove);

geo.fem=fem;
geo.boundnameAPmove = boundnameAPmove;
