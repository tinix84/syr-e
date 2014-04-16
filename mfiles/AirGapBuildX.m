%% 2013/08/21 MG contruction of arc and line of the airgap, bordo mobile dwaun apart...

function FemmProblem=AirGapBuildX(FemmProblem,Qs,ps,p,g,pc,xr,res_traf,groupStat,groupRot,boundnameAPg1,boundnameAPg2,boundnameAPg3)

[xArcStat1,yArcStat1] = rot_point(xr+2/3*g,0,-pc*pi/180);
[xArcStat2,yArcStat2] = rot_point(xr+2/3*g,0,(2*Qs-1)*pc*pi/180);

outernodes=[xArcStat1,yArcStat1;xArcStat2,yArcStat2];
[FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2));
[FemmProblem, ~] = addarcsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),...
                                                2*pc*Qs, 'MaxSegDegrees', res_traf,...
                                                'InGroup',groupStat);

% linee al traferro arco traf --> statore
[xTrafStat1,yTrafStat1]=rot_point(xr+g,0,-pc*pi/180);
[xTrafStat2,yTrafStat2]=rot_point(xr+g,0,(2*Qs-1)*pc*pi/180);
outernodes=[xTrafStat1,yTrafStat1;xTrafStat2,yTrafStat2];
[FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2));
[FemmProblem,id] = addsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),'InGroup',groupStat);
FemmProblem.Segments(id).BoundaryMarker = boundnameAPg3;

% mi_drawline(xArcStat1,yArcStat1,xTrafStat1,yTrafStat1);
% mi_selectsegment(mean([xArcStat1 xTrafStat1]),mean([yArcStat1 yTrafStat1]));
% mi_setsegmentprop('APg3', res_traf, 0, 0, groupStat);
% mi_selectnode(xTrafStat1,yTrafStat1); mi_setnodeprop('None',groupStat);
% mi_selectnode(xArcStat1,yArcStat1); mi_setnodeprop('None',groupStat);
% mi_clearselected;

outernodes=[xTrafStat1,yTrafStat1;xTrafStat2,yTrafStat2];
[FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2));
[FemmProblem,id] = addsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),'InGroup',groupStat);
FemmProblem.Segments(id).BoundaryMarker = boundnameAPg3;

% mi_drawline(xArcStat2,yArcStat2,xTrafStat2,yTrafStat2);
% mi_selectsegment(mean([xArcStat2 xTrafStat2]),mean([yArcStat2 yTrafStat2]));
% mi_setsegmentprop('APg3', res_traf, 0, 0, groupStat);
% mi_selectnode(xTrafStat2,yTrafStat2); mi_setnodeprop('None',groupStat);
% mi_selectnode(xArcStat2,yArcStat2); mi_setnodeprop('None',groupStat);
% mi_clearselected;

[xAirTrafSt,yAirTrafSt] = rot_point(xr+5/6*g,0,0.5*ps*180/p*pi/180);
FemmProblem = addblocklabel_mfemm(FemmProblem, xAirTrafSt,yAirTrafSt, ...
                                  'BlockType', 'Air', ...
                                  'MaxArea', res_traf,...
                                  'InGroup', groupStat);
                              
% mi_addblocklabel(xAirTrafSt,yAirTrafSt);
% mi_selectlabel(xAirTrafSt,yAirTrafSt);
% mi_setblockprop('Air', 0, res_traf, 'None', 0, groupStat, 1);
% mi_clearselected;
% linee al traferro arco traf --> stat bordo mobile

[xArcTraf1,yArcTraf1] = rot_point(xr+1/3*g,0,-pc*pi/180);
[xArcTraf2,yArcTraf2] = rot_point(xr+1/3*g,0,(2*Qs-1)*pc*pi/180);
outernodes=[xArcStat1,yArcStat1;xArcTraf1,yArcTraf1];
[FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2));
[FemmProblem,id] = addsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),'InGroup',groupStat);
FemmProblem.Segments(id).BoundaryMarker = boundnameAPg2;

% mi_drawline(xArcStat1,yArcStat1,xArcTraf1,yArcTraf1);
% mi_selectsegment(mean([xArcTraf1 xArcStat1]),mean([yArcTraf1 yArcStat1]));
% mi_setsegmentprop('APg2', res_traf, 0, 0, groupStat);
% mi_selectnode(xArcTraf1,yArcTraf1); mi_setnodeprop('None',groupStat);
% mi_clearselected;

outernodes=[xArcStat2,yArcStat2;xArcTraf2,yArcTraf2];
[FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2));
[FemmProblem,id] = addsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),'InGroup',groupStat);
FemmProblem.Segments(id).BoundaryMarker = boundnameAPg2;

% mi_drawline(xArcStat2,yArcStat2,xArcTraf2,yArcTraf2);
% mi_selectsegment(mean([xArcTraf2 xArcStat2]),mean([yArcTraf2 yArcStat2]));
% mi_setsegmentprop('APg2', res_traf, 0, 0, groupStat);
% mi_selectnode(xArcTraf2,yArcTraf2); mi_setnodeprop('None',groupStat);
% mi_clearselected;

[xAirTrafTr,yAirTrafTr] = rot_point(xr+3/6*g,0,0.5*ps*180/p*pi/180);
FemmProblem = addblocklabel_mfemm(FemmProblem, xAirTrafTr,yAirTrafTr, ...
                                  'BlockType', 'Air', ...
                                  'MaxArea', res_traf,...
                                  'InGroup', groupStat);
% mi_addblocklabel(xAirTrafTr,yAirTrafTr);
% mi_selectlabel(xAirTrafTr,yAirTrafTr);
% mi_setblockprop('Air', 0, res_traf, 'None', 0, groupStat, 1);
% mi_clearselected;
% disegno linee di rotore al --> traferro

[xArcRot1,yArcRot1] = rot_point(xr+1/3*g,0,0);
[xArcRot2,yArcRot2] = rot_point(xr+1/3*g,0,ps*180/p*pi/180);
[xRot1,yRot1] = rot_point(xr,0,0);
[xRot2,yRot2] = rot_point(xr,0,ps*180/p*pi/180);

outernodes=[xRot1,yRot1;xArcRot1,yArcRot1];
[FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2));
[FemmProblem,id] = addsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),'InGroup',groupRot);
FemmProblem.Segments(id).BoundaryMarker = boundnameAPg1;

% mi_drawline(xRot1,yRot1,xArcRot1,yArcRot1);
% mi_selectsegment(mean([xRot1 xArcRot1]),mean([yRot1 yArcRot1]));
% mi_setsegmentprop('APg1', res_traf, 0, 0, groupRot);
% mi_selectnode(xArcRot1,yArcRot1); mi_setnodeprop('None',groupRot);
% mi_clearselected;
outernodes=[xRot2,yRot2;xArcRot2,yArcRot2];
[FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2));
[FemmProblem,id] = addsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),'InGroup',groupRot);
FemmProblem.Segments(id).BoundaryMarker = boundnameAPg1;

% mi_drawline(xRot2,yRot2,xArcRot2,yArcRot2);
% mi_selectsegment(mean([xRot2 xArcRot2]),mean([yRot2 yArcRot2]));
% mi_setsegmentprop('APg1', res_traf, 0, 0, groupRot);
% mi_selectnode(xArcRot2,yArcRot2); mi_setnodeprop('None',groupRot);
% mi_clearselected;

[xAirTrafRot,yAirTrafRot] = rot_point(xr+1/6*g,0,0.5*ps*180/p*pi/180);
FemmProblem = addblocklabel_mfemm(FemmProblem, xAirTrafRot,yAirTrafRot, ...
                                  'BlockType', 'Air', ...
                                  'MaxArea', res_traf,...
                                  'InGroup', groupRot);
%                               
% mi_addblocklabel(xAirTrafRot,yAirTrafRot);
% mi_selectlabel(xAirTrafRot,yAirTrafRot);
% mi_setblockprop('Air', 0, res_traf, 'None', 0, groupRot, 1);
% mi_clearselected;

end