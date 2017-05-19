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
        
        if isfield(temp,'xPMso')
            xPMso = temp.xPMso;
            yPMso = temp.yPMso;
            xPMsi = temp.xPMsi;
            yPMsi = temp.yPMsi;
        end
        %% find the center of permanent magnet
        PMBaricentro = [];
        
        hybrid = geo.hybrid;
        if hybrid == 0
            %% regular or rounded shape
            if seg~=1
                % find the center of each PM segment; defined as the middle
                % point of diagonal;
                xPMcenter = (xPMi + xPMso(1))/2;
                yPMcenter = (yPMi + yPMso(1))/2;
                xPMcenter = [xPMcenter,(xPMsi(1:end-1) + xPMso(2:end))/2];
                yPMcenter = [yPMcenter,(yPMsi(1:end-1) + yPMso(2:end))/2];
                % duplicate to full pole
                if mod(seg,2)==1
                    %% odd segments
                    xPMcenter = [xPMcenter,geo.r-geo.lm*0.5];
                    yPMcenter = [yPMcenter,0];
                    xPMcenter = [xPMcenter,xPMcenter(1:end-1)];
                    yPMcenter = [yPMcenter,-yPMcenter(1:end-1)];
                else
                    %% even segments
                    xPMcenter = [xPMcenter,xPMcenter];
                    yPMcenter = [yPMcenter,-yPMcenter];
                end
                PMBaricentro = [PMBaricentro; xPMcenter',yPMcenter',codMatPM*ones(seg,1),res*ones(seg,1),ones(seg,1),zeros(seg,1),zeros(seg,1),zeros(seg,1)];
            else
                PMBaricentro = [mean([xPMci,xPMco]),0,codMatPM,res,1,0,0,0];
            end
        else
            %% hybrid shape
            PM_angle = geo.PM_angle;
            Fe_angle = geo.Fe_angle;
            [xFecenter1,yFecenter1] = rot_point(mean([xPMci,xPMco]),0,(PM_angle+Fe_angle/2)*pi/180);
            [xFecenter2,yFecenter2] = rot_point(mean([xPMci,xPMco]),0,(-PM_angle-Fe_angle/2)*pi/180);
            
            PMBaricentro = [mean([xPMci,xPMco]),0,codMatPM,res,1,0,0,0;
                xFecenter1,yFecenter1,codMatPM,res,1,0,0,0;
                xFecenter2,yFecenter2,codMatPM,res,1,0,0,0;];
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
        
        % Chao delete air zone
        Aircentro = [];
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
        
        if (0)  % put 1 to define the magnet direction parallel to q axis
            warning('magnet direction is parallel!!!')
            [r,c]=size(xmag);
            xmag=ones(r,c);
            ymag=zeros(r,c);
            zmag=zeros(r,c);
        end
        
        % barrier block centers
        if isempty(xc)
            BarCenter=[];
        else
            
            if(sum(geo.BarFillFac~=0)&&strcmp('Circular',geo.RotType))
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%just for circular type with BarFillFac~=0
                XpontRadBarSx=temp.XpontRadBarSx;
                YpontRadBarSx=temp.YpontRadBarSx;
                XpontRadBarDx=temp.XpontRadBarDx;
                YpontRadBarDx=temp.YpontRadBarDx;
                x0=geo.x0;
                Rbeta=sqrt((x0-xc).^2+yc.^2);
                Beta=atan((yc)./(x0-xc));
                if isfield(temp,'xxD1k')
                    Rbeta=x0-temp.Bx0;
                    Beta1=atan(temp.yyD1k./(x0-temp.xxD1k));
                    Beta2=atan(temp.yyD2k./(x0-temp.xxD2k));
                    Beta=zeros(1,length(Beta));
                    for ii=1:length(Beta)
                        if Beta1(ii)<Beta2(ii)
                            Beta(ii)=Beta1(ii);
                        else
                            Beta(ii)=Beta2(ii);
                        end
                    end
                end
                % Beta=atan(temp.Y5./(x0-temp.X5));
                Beta_min_dx=atan(temp.YpontRadDx./(x0-temp.XpontRadDx));
                Beta_min_sx=atan(temp.YpontRadSx./(x0-temp.XpontRadDx));
                Beta_min=zeros(1,length(Beta_min_dx));
                for ii=1:length(Beta_min_dx)
                    if Beta_min_dx(ii)<Beta_min_sx(ii)
                        Beta_min(ii)=Beta_min_sx(ii);
                    else
                        Beta_min(ii)=Beta_min_dx(ii);
                    end
                    if isnan(Beta_min(ii))
                        Beta_min(ii)=0;
                    end
                end
                Beta_half=(Beta+Beta_min)/2;
                Beta_mid=(Beta+Beta_min)/2;
                if isfield(temp,'xPMi')
                    xPMi=temp.xPMi;
                    yPMi=temp.yPMi;
                    xPMe=temp.xPMe;
                    yPMe=temp.yPMe;
                    Xmidbar=[];
                    Ymidbar=[];
                    for ii=1:length(xPMi)
                        Xmidbar=[Xmidbar,xPMi(ii),xPMe(ii)];
                        Ymidbar=[Ymidbar,yPMi(ii),yPMe(ii)];
                    end
                else
                    Xmidbar=(x0-Rbeta.*cos(Beta_mid));
                    Ymidbar=Rbeta.*sin(Beta_mid);
                end
                Xmidbar=[Xmidbar';Xmidbar']; Ymidbar=[Ymidbar';-Ymidbar'];
                if isfield(temp,'xpont')
                    xair=(temp.xpont+temp.xxD1k+temp.xxD2k)/3;
                    yair=(temp.ypont+temp.yyD1k+temp.yyD2k)/3;
                    xair=[xair';xair'];
                    yair=[yair';-yair'];
                else
                    yc_new=yc+([geo.B2k-geo.B1k]/4);
                    xair=[xc';xc'];
                    yair=[yc_new';-yc_new'];
                end
                %xmag=0*ones(1,length(xc)); ymag=0*ones(1,length(xc));zmag=0*ones(1,length(xc));
                
                xmag=[xmag';xmag']; ymag=[ymag';-ymag']; zmag=[zmag';zmag'];
                BarCenter=[];
                BarCenter=[BarCenter;Xmidbar,Ymidbar,codMatBar*ones(length(Xmidbar),1),res*ones(length(Xmidbar),1),1*ones(length(Xmidbar),1),xmag,ymag,zmag;...
                    xair,yair,codMatAir*ones(length(xair),1),res*ones(length(xair),1),1*ones(length(xair),1),zeros(length(xair),1),zeros(length(xair),1),zeros(length(xair),1)];
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
                Temp=[Temp;xtemp,ytemp,BarCenter(:,3),BarCenter(:,4),BarCenter(:,5),cos(magdir_tmp),sin(magdir_tmp),BarCenter(:,8)];
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
