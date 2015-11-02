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

function BarCenter = defineBlockCenters(temp,fem,geo)

% coordinates of center points of all FEMM blocks,
% where block labels will be eventually placed

res=fem.res;
ps = geo.ps;
p = geo.p;
Ar = geo.Ar;
nlay = geo.nlay;
phi = geo.phi/p;   % angle range of permanent magnet

% Material codes: refer to BLKLABELS.materials positions
codMatFe    = 5;
codMatBar   = 6;
codMatShaft = 7;
codMatPM    = 6;
codMatAir   = 1;

switch geo.RotType
    
    case 'SPM'
        
        xPMci = temp.xPMci;
        xPMco = temp.xPMco;
        % find the center of permanent magnet
        PMBaricentro = [mean([xPMci,xPMco]),0,codMatPM,res,1,0,0,0];
        % Replicate poles ps-1 times
        Temp=[]; kk=1;
        if not(isempty(PMBaricentro))
            while kk<=ps-1
                [xtemp,ytemp]=rot_point(PMBaricentro(:,1),PMBaricentro(:,2),kk*180/p*pi/180);
                Temp = [Temp; xtemp,ytemp,codMatPM,res,1,1,1,0];
                kk=kk+1;
            end
            PMBaricentro =[PMBaricentro;Temp];
            clear Temp xtemp ytemp;
        end
        %% find the center of air zone
        [xAir1,yAir1] = rot_point(mean([xPMci,xPMco]),0,(45/p+phi/4)*pi/180);
        [xAir2,yAir2] = rot_point(mean([xPMci,xPMco]),0,-(45/p+phi/4)*pi/180);
        Aircentro1 = [xAir1,yAir1,codMatAir,res,1,0,0,0];
        Aircentro2 = [xAir2,yAir2,codMatAir,res,1,0,0,0];
        Aircentro = [Aircentro1;Aircentro2];
        % Replicate poles ps-1 times
        Temp=[]; kk=1;
        if not(isempty(Aircentro))
            while kk<=ps-1
                [xtemp1,ytemp1]=rot_point(Aircentro1(:,1),Aircentro1(:,2),kk*180/p*pi/180);
                Temp = [Temp; xtemp1,ytemp1,codMatAir,res,1,0,0,0];
                [xtemp2,ytemp2]=rot_point(Aircentro2(:,1),Aircentro2(:,2),kk*180/p*pi/180);
                Temp = [Temp; xtemp2,ytemp2,codMatAir,res,1,0,0,0];
                kk=kk+1;
            end
            Aircentro =[Aircentro;Temp];
            clear Temp xtemp ytemp;
        end
        %% find the center of iron
        RotBaricentro=[mean([Ar,xPMci]),0,codMatFe,res,1,NaN,NaN,NaN];
        %% find the center of shaft
        ShaftBaricentro=[mean([0,Ar]),0,codMatShaft,res,1,NaN,NaN,NaN];
        
        BarCenter=[PMBaricentro;Aircentro;RotBaricentro;ShaftBaricentro];
        
    otherwise
        
        xc = temp.xc;
        yc = temp.yc;
        
        % Barriers:
        % magnetization direction
        xmag=temp.xmag; ymag=temp.ymag; zmag=temp.zmag;
        
        % barrier block centers
        if isempty(xc)
            BarCenter=[];
        else
            xc=[xc';xc']; yc=[yc';-yc'];
            xmag=[xmag';xmag']; ymag=[ymag';-ymag']; zmag=[zmag';zmag'];
            BarCenter=[xc,yc,codMatBar*ones(length(xc),1),res*ones(length(xc),1),1*ones(length(xc),1),xmag,ymag,zmag];
            magdir=atan2(ymag,xmag);
        end
        
        % Replicate poles ps-1 times
        Temp=[]; kk=1;
        if not(isempty(BarCenter))
            while kk<=ps-1
                [xtemp,ytemp]=rot_point(BarCenter(:,1),BarCenter(:,2),kk*180/p*pi/180);
                magdir_tmp=magdir+(kk*pi/p-eps)+(cos((kk-1)*pi)+1)/2*pi;
                Temp=[Temp;xtemp,ytemp,codMatBar*ones(length(xc),1),res*ones(length(xc),1),1*ones(length(xc),1),cos(magdir_tmp),sin(magdir_tmp),zmag];
                kk=kk+1;
            end
            BarCenter=[BarCenter;Temp];
            clear Temp xtemp ytemp;
        end
                
        % rotor iron
        XBan1sx=geo.B1k;
        RotBaricentro=[mean([Ar,XBan1sx(nlay)]),0,codMatFe,res,1,NaN,NaN,NaN];
        
        % shaft
        ShaftBaricentro=[mean([0,Ar]),0,codMatShaft,res,1,NaN,NaN,NaN];
        
        BarCenter=[BarCenter;RotBaricentro;ShaftBaricentro];
        
end
