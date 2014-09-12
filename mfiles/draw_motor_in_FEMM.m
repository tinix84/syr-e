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

function [geo] = draw_motor_in_FEMM(geo,eval_type)

% draw_motor_in_FEMM.m
% builds the .fem motor model with the rotor in position zero

% input: see function argument
% output:
% - updated geo (individual for the machine)
% - mot0.fem (th_m = 0, i123 = [0,0,0])

fem = dimMesh(geo,eval_type);
% fem.res_traf=0.5;
% fem.res=6;

filename = 'mot0.fem';
opendocument('empty_case.fem');
mi_probdef(0,'millimeters','planar',1e-8,geo.l,25);

% calc winding factor (kavv) and rotor offset (phase1_offset)
[kavv, phase1_offset] = calcola_kavv_th0(geo.avv,geo.ns*geo.p,size(geo.avv,1),geo.p);

% offset angle for coordinate transformations
th_m0 = 0;                              % rotor position [mec deg]
geo.th0 = th_m0*geo.p - phase1_offset;  % d- to alpha-axis offset [elt deg]

% machine periodicity (t) and number of poles to be simulated (ps)
Q=geo.ns*geo.p;
t=gcd(geo.ns*geo.p,geo.p);  % periodicity
if ((6*t/Q)>1)
    ps=2*geo.p/t;
    Qs=Q/t;
else
    ps=geo.p/t;
    Qs=Q/2/t;
end

geo.Qs=Qs;  % # of simulated slots
geo.ps=ps;  % # of simulated poles

% Boundary conditions
if (rem(geo.ps,2)==0)
    periodicity=4;
else
    periodicity=5;
end

%% definition of the boundary conditions
% inner and outer circles
mi_addboundprop('A=0', 0, 0, 0, 0, 0, 0, 0, 0, 0);
% Periodicity or Anti-Periodicity (2 x rotor + 2 x stator + 3 x airgap + 1 x APmove) APmove is the sliding contour
mi_addboundprop('APr1', 0, 0, 0, 0, 0, 0, 0, 0, periodicity);
mi_addboundprop('APr2', 0, 0, 0, 0, 0, 0, 0, 0, periodicity);
mi_addboundprop('APg1', 0, 0, 0, 0, 0, 0, 0, 0, periodicity);
mi_addboundprop('APg2', 0, 0, 0, 0, 0, 0, 0, 0, periodicity);
mi_addboundprop('APg3', 0, 0, 0, 0, 0, 0, 0, 0, periodicity);
mi_addboundprop('APmove', 0, 0, 0, 0, 0, 0, 0, 0, periodicity);
mi_addboundprop('APs1', 0, 0, 0, 0, 0, 0, 0, 0, periodicity);
mi_addboundprop('APs2', 0, 0, 0, 0, 0, 0, 0, 0, periodicity);

%% ROTOR
% build the matrices which describe the rotor
geo.x0 = geo.xr/cos(pi/2/geo.p);
ROTmatr;
BLKLABELS.rotore=BLKLABELSrot;

% draw lines and arcs
draw_lines_arches(rotore2,2,fem.res);

% assign block labels
assign_block_prop_rot(BLKLABELS,geo,fem,2);

% prevents flux barrier overlapping in high speed ISeg machines (with radial ribs)
% ISeg_HS_bug_fix;

% assign boundary conditions
for ii=1:2
    mi_selectsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2));
    if (BLKLABELSrot.boundary(ii,3)==10)
        mi_setsegmentprop('APr1', 0, 1, 0, 2);
        mi_clearselected;
    elseif(BLKLABELSrot.boundary(ii,3)==0)
        mi_selectarcsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2))
        mi_setarcsegmentprop(fem.res, 'A=0', 0, 2);
        mi_clearselected;
    end
end
for ii=3:4   
    mi_selectsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2));
    if (BLKLABELSrot.boundary(ii,3)==10)
        mi_setsegmentprop('APr2', 0, 1, 0, 2);
        mi_clearselected;
    elseif(BLKLABELSrot.boundary(ii,3)==0)
        mi_selectarcsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2))
        mi_setarcsegmentprop(fem.res, 'A=0', 0, 2);
        mi_clearselected;
    end
end

%% STATOR
% builds the matrixes which describe the stator
STATmatr;
BLKLABELS.statore=BLKLABELSstat;

% draw lines and arcs
draw_lines_arches(statore,1,fem.res);

% assign block labels
assign_block_prop_stat(BLKLABELS,geo,fem,1) % assegna materiali;

% assign boundary conditions
BLKLABELSstat=BLKLABELS.statore;
for ii=1:size(BLKLABELSstat.boundary,1)    
    mi_selectsegment(BLKLABELSstat.boundary(ii,1),BLKLABELSstat.boundary(ii,2));
    if (BLKLABELSstat.boundary(ii,3)==10)
        mi_setsegmentprop('APs1', 0, 1, 0, 1);
        mi_clearselected;
    elseif(BLKLABELSstat.boundary(ii,3)==0)
        mi_selectarcsegment(BLKLABELSstat.boundary(ii,1),BLKLABELSstat.boundary(ii,2));
        mi_setarcsegmentprop(fem.res, 'A=0', 0, 1);
        mi_clearselected;
    end
end

%% airgap (group 20)
    AirGapBuild(Qs,ps,geo.p,geo.g,360/(ns*geo.p)/2,geo.xr,fem.res_traf,1,2);
    draw_airgap_arc_with_mesh(geo,th_m0,fem);

geo.fem=fem;
mi_saveas(filename); % saves the file with name ’filename’.
