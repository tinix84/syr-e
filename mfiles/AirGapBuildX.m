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

%% 2013/08/21 MG contruction of arc and line of the airgap, bordo mobile dwaun apart...

function FemmProblem=AirGapBuildX(FemmProblem,Qs,ps,p,g,pc,xr,res_traf,groupStat,groupRot,boundnameAPg1,boundnameAPg2,boundnameAPg3)

[xArcStat1,yArcStat1] = rot_point(xr+2/3*g,0,-pc*pi/180);
[xArcStat2,yArcStat2] = rot_point(xr+2/3*g,0,(2*Qs-1)*pc*pi/180);

outernodes=[xArcStat1,yArcStat1;xArcStat2,yArcStat2];
[FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2),'InGroup',groupStat);
[FemmProblem, ~] = addarcsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),...
                                                2*pc*Qs, 'MaxSegDegrees', res_traf,...
                                                'InGroup',groupStat);

% linee al traferro arco traf --> statore
[xTrafStat1,yTrafStat1]=rot_point(xr+g,0,-pc*pi/180);
[xTrafStat2,yTrafStat2]=rot_point(xr+g,0,(2*Qs-1)*pc*pi/180);
outernodes=[xArcStat1,yArcStat1;xTrafStat1,yTrafStat1];
[FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2),'InGroup',groupStat);
[FemmProblem,id] = addsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),'InGroup',groupStat);
FemmProblem.Segments(id).BoundaryMarker = boundnameAPg3;

outernodes=[xArcStat2,yArcStat2;xTrafStat2,yTrafStat2];
[FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2),'InGroup',groupStat);
[FemmProblem,id] = addsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),'InGroup',groupStat);
FemmProblem.Segments(id).BoundaryMarker = boundnameAPg3;

[xAirTrafSt,yAirTrafSt] = rot_point(xr+5/6*g,0,0.5*ps*180/p*pi/180);
FemmProblem = addblocklabel_mfemm(FemmProblem, xAirTrafSt,yAirTrafSt, ...
                                  'BlockType', 'Air', ...
                                  'MaxArea', res_traf,...
                                  'InGroup', groupStat);
                              
[xArcTraf1,yArcTraf1] = rot_point(xr+1/3*g,0,-pc*pi/180);
[xArcTraf2,yArcTraf2] = rot_point(xr+1/3*g,0,(2*Qs-1)*pc*pi/180);
outernodes=[xArcStat1,yArcStat1;xArcTraf1,yArcTraf1];
[FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2),'InGroup',groupStat);
[FemmProblem,id] = addsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),'InGroup',groupStat);
FemmProblem.Segments(id).BoundaryMarker = boundnameAPg2;

outernodes=[xArcStat2,yArcStat2;xArcTraf2,yArcTraf2];
[FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2),'InGroup',groupStat);
[FemmProblem,id] = addsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),'InGroup',groupStat);
FemmProblem.Segments(id).BoundaryMarker = boundnameAPg2;


[xAirTrafTr,yAirTrafTr] = rot_point(xr+3/6*g,0,0.5*ps*180/p*pi/180);
FemmProblem = addblocklabel_mfemm(FemmProblem, xAirTrafTr,yAirTrafTr, ...
                                  'BlockType', 'Air', ...
                                  'MaxArea', res_traf,...
                                  'InGroup', groupStat);

% disegno linee di rotore al --> traferro

[xArcRot1,yArcRot1] = rot_point(xr+1/3*g,0,0);
[xArcRot2,yArcRot2] = rot_point(xr+1/3*g,0,ps*180/p*pi/180);
[xRot1,yRot1] = rot_point(xr,0,0);
[xRot2,yRot2] = rot_point(xr,0,ps*180/p*pi/180);

outernodes=[xRot1,yRot1;xArcRot1,yArcRot1];
[FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2),'InGroup',groupRot);
[FemmProblem,id] = addsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),'InGroup',groupRot);
FemmProblem.Segments(id).BoundaryMarker = boundnameAPg1;

outernodes=[xRot2,yRot2;xArcRot2,yArcRot2];
[FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2),'InGroup',groupRot);
[FemmProblem,id] = addsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),'InGroup',groupRot);
FemmProblem.Segments(id).BoundaryMarker = boundnameAPg1;

[xAirTrafRot,yAirTrafRot] = rot_point(xr+1/6*g,0,0.5*ps*180/p*pi/180);
FemmProblem = addblocklabel_mfemm(FemmProblem, xAirTrafRot,yAirTrafRot, ...
                                  'BlockType', 'Air', ...
                                  'MaxArea', res_traf,...
                                  'InGroup', groupRot);

end