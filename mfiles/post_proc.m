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


l = geo.l;

% load phase flux linkages
for ii=0:(n3phase-1) %AS
    temp_out = mo_getcircuitproperties(phase_name{3*ii+1});
    temp_out = temp_out - mo_getcircuitproperties(phase_name_neg{3*ii+1});
    f(3*ii+1) = temp_out(3) * 2 * p/ps; %ps number of poles in FEMM
    temp_out = mo_getcircuitproperties(phase_name{3*ii+2});
    temp_out = temp_out - mo_getcircuitproperties(phase_name_neg{3*ii+2});
    f(3*ii+2) = temp_out(3) * 2 * p/ps;
    temp_out = mo_getcircuitproperties(phase_name{3*ii+3});
    temp_out = temp_out - mo_getcircuitproperties(phase_name_neg{3*ii+3});
    f(3*ii+3) = temp_out(3) * 2 * p/ps;
end

% evaluate torque
% % T1 - from the innermost integration line (rot + gap/6)
% x = r + gap*1/6;
% if gradi_da_sim<=180
%     ang0 = th_m;
%     ang1 = gradi_da_sim+th_m;
%     [x1,y1] = rot_point(x,0,ang0*pi/180);
%     [x2,y2] = rot_point(x,0,ang1*pi/180);
%     mo_addcontour(x1,y1);
%     mo_addcontour(x2,y2);
%     mo_bendcontour(gradi_da_sim,0.5);
% else
%     ang0 = th_m;
%     ang1 = gradi_da_sim/2+th_m;
%     ang2 = gradi_da_sim+th_m;
%     [x1,y1] = rot_point(x,0,ang0*pi/180);
%     [x2,y2] = rot_point(x,0,ang1*pi/180);
%     [x3,y3] = rot_point(x,0,ang2*pi/180);
%     mo_addcontour(x1,y1);
%     mo_addcontour(x2,y2);
%     mo_bendcontour(gradi_da_sim/2,0.5);
%     mo_addcontour(x3,y3);
%     mo_bendcontour(gradi_da_sim/2,0.5);
% end

% T1 = mo_lineintegral(4);
% T1 = T1(1) * 2 * p/ps;
% mo_clearcontour();

% % T2 - from the outermost integration line (stat - gap/6)
% x = r + gap*5/6;
% if gradi_da_sim<=180
%     ang0 = -pc;
%     ang1 = gradi_da_sim-pc;
%     [x1,y1] = rot_point(x,0,ang0*pi/180);
%     [x2,y2] = rot_point(x,0,ang1*pi/180);
%     mo_addcontour(x1,y1);
%     mo_addcontour(x2,y2);
%     mo_bendcontour(gradi_da_sim,0.5);
% else
%     ang0 = -pc;
%     ang1 = gradi_da_sim/2-pc;
%     ang2 = gradi_da_sim-pc;
%     [x1,y1] = rot_point(x,0,ang0*pi/180);
%     [x2,y2] = rot_point(x,0,ang1*pi/180);
%     [x3,y3] = rot_point(x,0,ang2*pi/180);
%     mo_addcontour(x1,y1);
%     mo_addcontour(x2,y2);
%     mo_bendcontour(gradi_da_sim/2,0.5);
%     mo_addcontour(x3,y3);
%     mo_bendcontour(gradi_da_sim/2,0.5);
% end
% T2 = mo_lineintegral(4);
% T2 = T2(1) * 2 * p/ps;
% mo_clearcontour();

% % T3 - from an intermediate line (rot + gap/2)
% x = r + gap*1/2;
% if gradi_da_sim<=180
%     ang0 = -pc;
%     ang1 = gradi_da_sim-pc;
%     [x1,y1] = rot_point(x,0,ang0*pi/180);
%     [x2,y2] = rot_point(x,0,ang1*pi/180);
%     mo_addcontour(x1,y1);
%     mo_addcontour(x2,y2);
%     mo_bendcontour(gradi_da_sim,0.5);
% else
%     ang0 = -pc;
%     ang1 = gradi_da_sim/2-pc;
%     ang2 = gradi_da_sim-pc;
%     [x1,y1] = rot_point(x,0,ang0*pi/180);
%     [x2,y2] = rot_point(x,0,ang1*pi/180);
%     [x3,y3] = rot_point(x,0,ang2*pi/180);
%     mo_addcontour(x1,y1);
%     mo_addcontour(x2,y2);
%     mo_bendcontour(gradi_da_sim/2,0.5);
%     mo_addcontour(x3,y3);
%     mo_bendcontour(gradi_da_sim/2,0.5);
% end
% T3 = mo_lineintegral(4);
% T3 = T3(1) * 2 * p/ps;
% mo_clearcontour();

% dq flux linkaged
for ik=0:(n3phase-1) %AS
    fdq = abc2dq(f(3*ik+1),f(3*ik+2),f(3*ik+3),th(jj)*pi/180,n3phase,ik);
%     flusso(ik+1,jj)=fdq(1)+j*fdq(2);
    fd_temp(ik+1,jj)=fdq(1);
    fq_temp(ik+1,jj)=fdq(2);
end

fd=mean(fd_temp(:,jj));
fq=mean(fq_temp(:,jj));

% block evaluation of torque and force
for ii=1:length(geo.BLKLABELS.rotore.xy(:,1))
    xB=geo.BLKLABELS.rotore.xy(ii,1);
    yB=geo.BLKLABELS.rotore.xy(ii,2);
    %[xB,yB]=rot_point(xB,yB,(th(jj) - th(1))/p*pi/180);
    [xB,yB]=rot_point(xB,yB,th_m*pi/180);
    mo_selectblock(xB,yB);
end

Tblock=mo_blockintegral(22)*2*p/ps;
%Fx=mo_blockintegral(18);
%Fy=mo_blockintegral(19);
mo_clearblock;

% Calcolo del Volume dei magneti VolPM - rev.Gallo 14/03/2018
if jj==1 %viene calcolato solo alla prima simulazione (rotore in posizione di partenza)
    flagPM=0;
    for ii=1:length(geo.BLKLABELS.rotore.xy(:,1))
        if geo.BLKLABELS.rotore.xy(ii,3)== 6 %lettura codice del materiale per individuare regioni di PM presenti nella struttura di rotore
            xA=geo.BLKLABELS.rotore.xy(ii,1);
            yA=geo.BLKLABELS.rotore.xy(ii,2);
            [xA,yA]=rot_point(xA,yA,th_m*pi/180); %rotazione delle coordinate di random position offset
            mo_selectblock(xA,yA);
            flagPM=1;
        end
    end
    if flagPM %flag per capire se area del magnete rettangolare è presente o no 
        VolPM=(2*geo.p*mo_blockintegral(10))/geo.ps; %Calcolo Volume magnete totale nel rotore [m3]
    else
        VolPM=0;
    end
    
    mo_clearblock();
end

% ang0 = th_m;
% ang1 = gradi_da_sim+th_m;
% angC = mean([ang0 ang1]);
% Frt=(Fx+j*Fy)*exp(-j*angC);
% string of SOL
% sol = [th(jj) id iq fdq(1) fdq(2) mean([T1,T2,T3])];
% %sol = [th(jj) id iq fdq(1) fdq(2) mean([T1,T2,T3]) T1 T2 T3];
% sol = [th(jj) id iq fdq(1) fdq(2) mean([T1,T2,T3]) mean([F1,F2,F3]) F1 F2 F3];




