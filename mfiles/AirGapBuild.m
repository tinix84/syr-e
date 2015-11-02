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

%% 2013/08/21 MG contruction of arc and line of the airgap, bordo mobile dwaun apart...

function AirGapBuild(Qs,ps,p,g,pc,r,res_traf,groupStat,groupRot,lm,BLKLABELSrot,RotType)

[xArcStat1,yArcStat1] = rot_point(r+2/3*g,0,-pc*pi/180);
[xArcStat2,yArcStat2] = rot_point(r+2/3*g,0,(2*Qs-1)*pc*pi/180);

mi_drawarc(xArcStat1,yArcStat1,xArcStat2,yArcStat2,2*pc*Qs,res_traf);
mi_selectarcsegment(xArcStat1,yArcStat1);
mi_setarcsegmentprop(res_traf, 'None', 0, groupStat);
mi_selectnode(xArcStat1,yArcStat1); mi_setnodeprop('None',groupStat);
mi_selectnode(xArcStat2,yArcStat2); mi_setnodeprop('None',groupStat);
mi_clearselected

% linee al traferro arco traf --> statore
[xTrafStat1,yTrafStat1]=rot_point(r+g,0,-pc*pi/180);
[xTrafStat2,yTrafStat2]=rot_point(r+g,0,(2*Qs-1)*pc*pi/180);
mi_drawline(xArcStat1,yArcStat1,xTrafStat1,yTrafStat1);
mi_selectsegment(mean([xArcStat1 xTrafStat1]),mean([yArcStat1 yTrafStat1]));
mi_setsegmentprop('APg3', res_traf, 0, 0, groupStat);
mi_selectnode(xTrafStat1,yTrafStat1); mi_setnodeprop('None',groupStat);
mi_selectnode(xArcStat1,yArcStat1); mi_setnodeprop('None',groupStat);
mi_clearselected;

mi_drawline(xArcStat2,yArcStat2,xTrafStat2,yTrafStat2);
mi_selectsegment(mean([xArcStat2 xTrafStat2]),mean([yArcStat2 yTrafStat2]));
mi_setsegmentprop('APg3', res_traf, 0, 0, groupStat);
mi_selectnode(xTrafStat2,yTrafStat2); mi_setnodeprop('None',groupStat);
mi_selectnode(xArcStat2,yArcStat2); mi_setnodeprop('None',groupStat);
mi_clearselected;

[xAirTrafSt,yAirTrafSt] = rot_point(r+5/6*g,0,0.5*ps*180/p*pi/180);
mi_addblocklabel(xAirTrafSt,yAirTrafSt);
mi_selectlabel(xAirTrafSt,yAirTrafSt);
mi_setblockprop('Air', 0, res_traf, 'None', 0, groupStat, 1);
mi_clearselected;
% linee al traferro arco traf --> stat bordo mobile
[xArcTraf1,yArcTraf1] = rot_point(r+1/3*g,0,-pc*pi/180);
[xArcTraf2,yArcTraf2] = rot_point(r+1/3*g,0,(2*Qs-1)*pc*pi/180);
mi_drawline(xArcStat1,yArcStat1,xArcTraf1,yArcTraf1);
mi_selectsegment(mean([xArcTraf1 xArcStat1]),mean([yArcTraf1 yArcStat1]));
mi_setsegmentprop('APg2', res_traf, 0, 0, groupStat);
mi_selectnode(xArcTraf1,yArcTraf1); mi_setnodeprop('None',groupStat);
mi_clearselected;

mi_drawline(xArcStat2,yArcStat2,xArcTraf2,yArcTraf2);
mi_selectsegment(mean([xArcTraf2 xArcStat2]),mean([yArcTraf2 yArcStat2]));
mi_setsegmentprop('APg2', res_traf, 0, 0, groupStat);
mi_selectnode(xArcTraf2,yArcTraf2); mi_setnodeprop('None',groupStat);
mi_clearselected;

[xAirTrafTr,yAirTrafTr] = rot_point(r+3/6*g,0,0.5*ps*180/p*pi/180);
mi_addblocklabel(xAirTrafTr,yAirTrafTr);
mi_selectlabel(xAirTrafTr,yAirTrafTr);
mi_setblockprop('Air', 0, res_traf, 'None', 0, groupStat, 1);
mi_clearselected;
% disegno linee di rotore al --> traferro

[xArcRot1,yArcRot1] = rot_point(r+1/3*g,0,0);
[xArcRot2,yArcRot2] = rot_point(r+1/3*g,0,ps*180/p*pi/180);
if strcmp(RotType,'SPM')
    [xRot1,yRot1] = rot_point(r-lm,0,0);
    [xRot2,yRot2] = rot_point(r-lm,0,ps*180/p*pi/180);
else
    [xRot1,yRot1] = rot_point(r,0,0);
    [xRot2,yRot2] = rot_point(r,0,ps*180/p*pi/180);
end
mi_clearselected;

mi_drawline(xRot1,yRot1,xArcRot1,yArcRot1);
mi_selectsegment(mean([xRot1 xArcRot1]),mean([yRot1 yArcRot1]));
mi_setsegmentprop('APg1', res_traf, 0, 0, groupRot);
mi_selectnode(xArcRot1,yArcRot1); mi_setnodeprop('None',groupRot);
mi_clearselected;

mi_drawline(xRot2,yRot2,xArcRot2,yArcRot2);
mi_selectsegment(mean([xRot2 xArcRot2]),mean([yRot2 yArcRot2]));
mi_setsegmentprop('APg1', res_traf, 0, 0, groupRot);
mi_selectnode(xArcRot2,yArcRot2); mi_setnodeprop('None',groupRot);
mi_clearselected;

[xAirTrafRot,yAirTrafRot] = rot_point(r+1/6*g,0,0.5*ps*180/p*pi/180);
mi_addblocklabel(xAirTrafRot,yAirTrafRot);
mi_selectlabel(xAirTrafRot,yAirTrafRot);
mi_setblockprop('Air', 0, res_traf, 'None', 0, groupRot, 1);
mi_clearselected;

end