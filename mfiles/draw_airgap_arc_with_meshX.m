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

% rebuilds the three arcs on the sliding surface between the stator and the
% rotor
% re-assigns the AP boundary conditions
% re-assigns the angular resolution for having a regular airgap mesh

function FemmProblem=draw_airgap_arc_with_meshX(FemmProblem,geo,th_m,fem,boundnameAPmove)

p  = geo.p;
xr = geo.xr;
g  = geo.g;
ns = geo.ns;
ps=geo.ps;
% Mezzo passo cava
pc = 360/(ns*p)/2;

res=fem.res_traf;

% si fa riferimento ai gradi meccanici della porzione di rot simulata
gradi_da_sim = 180/p*ps;

group = 20;

x0 = xr+1/3*g;
y0 = 0;

res_deg=res;
angoli_bordo_mobile = [-pc, th_m, (-pc + gradi_da_sim), (gradi_da_sim+th_m)];
angoli_bordo_mobile = sort(angoli_bordo_mobile);
delta_angoli_bordo_mobile = diff(angoli_bordo_mobile);

[x1,y1] = rot_point(x0,y0,angoli_bordo_mobile(1)*pi/180);
[x2,y2] = rot_point(x0,y0,angoli_bordo_mobile(2)*pi/180);
[x3,y3] = rot_point(x0,y0,angoli_bordo_mobile(3)*pi/180);
[x4,y4] = rot_point(x0,y0,angoli_bordo_mobile(4)*pi/180);
outernodes=[x1,y1;x2,y2];
[FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2),'InGroup',group);
[FemmProblem, id] = addarcsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),...
                                                delta_angoli_bordo_mobile(1),...
                                                'MaxSegDegrees', res,...
                                                'InGroup',group);
FemmProblem.ArcSegments(id).BoundaryMarker = boundnameAPmove;
% num_seg=delta_angoli_bordo_mobile(1)*pi/180*x0*res;
% mi_addarc(x1,y1,x2,y2,delta_angoli_bordo_mobile(1),res);
% [x,y] = rot_point(x1,y1,0.5*delta_angoli_bordo_mobile(1)*pi/180);
% mi_selectarcsegment(x,y);
% mi_setarcsegmentprop(res_deg, 'APmove', 0, group);

% num_seg=delta_angoli_bordo_mobile(2)*pi/180*x0*res/ps;
outernodes=[x2,y2;x3,y3];
[FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2),'InGroup',group);
[FemmProblem, id] = addarcsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),...
                                                delta_angoli_bordo_mobile(2),...
                                                'MaxSegDegrees', res,...
                                                'InGroup',group);

% num_seg=delta_angoli_bordo_mobile(3)*pi/180*x0/res;

outernodes=[x3,y3;x4,y4];
[FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2),'InGroup',group);
[FemmProblem, id] = addarcsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),...
                                                delta_angoli_bordo_mobile(3),...
                                                'MaxSegDegrees', res,...
                                                'InGroup',group);
FemmProblem.ArcSegments(id).BoundaryMarker = boundnameAPmove;

% mi_addarc(x3,y3,x4,y4,delta_angoli_bordo_mobile(3),res);
% [x,y] = rot_point(x3,y3,0.5*delta_angoli_bordo_mobile(3)*pi/180);
% mi_selectarcsegment(x,y);
% mi_setarcsegmentprop(res_deg, 'APmove', 0, group);
% mi_clearselected
end



