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
if strcmp(geo.RotType,'SPM')
    seg = geo.dx;                     % number of segments of magnet for SPM
end

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
        xPMo = temp.xPMo;
        xPMi = temp.xPMi;
        yPMo = temp.yPMo;
        yPMi = temp.yPMi;
        %% find the center of permanent magnet
        PMBaricentro = [];
        
        if seg~=1
            xPMcenter0 = mean([xPMo,xPMi]);
            yPMcenter0 = mean([yPMo,yPMi]);
            for jj = 1:seg
                [xPMcenter(jj),yPMcenter(jj)] = rot_point(xPMcenter0,yPMcenter0,-(phi/2*pi/180/seg+phi*pi/180/seg*(jj-1)));
                PMBaricentro = [PMBaricentro; xPMcenter(jj),yPMcenter(jj),codMatPM,res,1,0,0,0];
            end
        else
            PMBaricentro = [mean([xPMci,xPMco]),0,codMatPM,res,1,0,0,0];
        end
        % Replicate poles ps-1 times
        Temp=[]; kk=1;
        
        if not(isempty(PMBaricentro))
            while kk<=ps-1
                for ii = 1:seg
                    [xtemp(ii),ytemp(ii)]=rot_point(PMBaricentro(ii,1),PMBaricentro(ii,2),kk*180/p*pi/180);
                    Temp = [Temp; xtemp(ii),ytemp(ii),codMatPM,res,1,0,0,0];
                end
                kk=kk+1;
            end
            PMBaricentro =[PMBaricentro;Temp];
            clear Temp xtemp ytemp;
        end
        %% find the center of air zone
        [xAir1,yAir1] = rot_point(mean([xPMci,xPMco]),0,(45/p+phi/4)*pi/180);
        [xAir2,yAir2] = rot_point(mean([xPMci,xPMco]),0,-(45/p+phi/4)*pi/180);
        Aircentro1 = [xAir1,yAir1,codMatAir,res,1,NaN,NaN,0];
        Aircentro2 = [xAir2,yAir2,codMatAir,res,1,NaN,NaN,0];
        Aircentro = [Aircentro1;Aircentro2];
        % Replicate poles ps-1 times
        Temp=[]; kk=1;
        if not(isempty(Aircentro))
            while kk<=ps-1
                [xtemp1,ytemp1]=rot_point(Aircentro1(:,1),Aircentro1(:,2),kk*180/p*pi/180);
                Temp = [Temp; xtemp1,ytemp1,codMatAir,res,1,NaN,NaN,0];
                [xtemp2,ytemp2]=rot_point(Aircentro2(:,1),Aircentro2(:,2),kk*180/p*pi/180);
                Temp = [Temp; xtemp2,ytemp2,codMatAir,res,1,NaN,NaN,0];
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
        
        % Barriers: magnetization direction
        xmag=temp.xmag; ymag=temp.ymag; zmag=temp.zmag;
        
        % barrier block centers
        if isempty(xc)
            BarCenter=[];
        else
            
%             a=(geo.BarFillFac~=0)&&strcmp('Circular',geo.RotType);
%             if (a==0)
            if(sum(geo.BarFillFac~=0)&&strcmp('Circular',geo.RotType))
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%just for circular type with BarFillFac~=0
                XpontRadBarSx=temp.XpontRadBarSx;
                YpontRadBarSx=temp.YpontRadBarSx;
                XpontRadBarDx=temp.XpontRadBarDx;
                YpontRadBarDx=temp.YpontRadBarDx;
                x0=geo.x0;
                Rbeta=sqrt((x0-xc).^2+yc.^2);
                Beta=atan((yc)./(x0-xc));
                Beta_mid=(Beta.*geo.BarFillFac)/2;            %%Determine center of magnet area
                Xmidbar=(x0-Rbeta.*cos(Beta_mid));
                Ymidbar=Rbeta.*sin(Beta_mid);
                
                Xmidbar=[Xmidbar';Xmidbar']; Ymidbar=[Ymidbar';-Ymidbar'];
                xc=xc;yc_new=yc+([geo.B2k-geo.B1k]/4);
                %xmag=0*ones(1,length(xc)); ymag=0*ones(1,length(xc));zmag=0*ones(1,length(xc));
                xc=[xc';xc']; yc_new=[yc_new';-yc_new'];
                xmag=[xmag';xmag']; ymag=[ymag';-ymag']; zmag=[zmag';zmag'];
                BarCenter=[];
                BarCenter=[BarCenter;Xmidbar,Ymidbar,codMatBar*ones(length(xc),1),res*ones(length(xc),1),1*ones(length(xc),1),xmag,ymag,zmag;...
                    xc,yc_new,codMatAir*ones(length(xc),1),res*ones(length(xc),1),1*ones(length(xc),1),xmag,ymag,zmag];
                % magdir=atan2(ymag,xmag);
                magdir=atan2(BarCenter(:,7),BarCenter(:,6)); % it means : magdir=atan2(ymag,xmag);
                % keyboard
            else
                xc=[xc';xc']; yc=[yc';-yc'];
                xmag=[xmag';xmag']; ymag=[ymag';-ymag']; zmag=[zmag';zmag'];
                BarCenter=[xc,yc,codMatBar*ones(length(xc),1),res*ones(length(xc),1),1*ones(length(xc),1),xmag,ymag,zmag];
                magdir=atan2(ymag,xmag);
                
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Replicate poles ps-1 times
        Temp=[]; kk=1;
        if not(isempty(BarCenter))
            while kk<=ps-1
                [xtemp,ytemp]=rot_point(BarCenter(:,1),BarCenter(:,2),kk*180/p*pi/180);
                magdir_tmp=magdir+(kk*pi/p-eps)+(cos((kk-1)*pi)+1)/2*pi;
                Temp=[Temp;xtemp,ytemp,codMatBar*ones(length(BarCenter(:,1)),1),res*ones(length(BarCenter(:,1)),1),1*ones(length(BarCenter(:,1)),1),cos(magdir_tmp),sin(magdir_tmp),BarCenter(:,8)];
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
