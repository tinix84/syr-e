
% Copyright 2018
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

function [geo,mat,temp] = MagnetFullFill_Seg2_fun(geo,mat,temp)

% assign the condition for flux barrier division

nlay = geo.nlay;

xxD1k = temp.xxD1k;
yyD1k = temp.yyD1k;
xxD2k = temp.xxD2k;
yyD2k = temp.yyD2k;
XpBar1 = temp.XpBar1;
YpBar1 = temp.YpBar1;
XpBar2 = temp.XpBar2;
YpBar2 = temp.YpBar2;
XpMag1B1 = temp.XpMag1B1;
YpMag1B1 = temp.YpMag1B1;
xpont = temp.xpont;
ypont = temp.ypont;
B1k = temp.B1k;
B2k = temp.B2k;
YpontRadBarDx = temp.YpontRadBarDx;
YpontRadBarSx = temp.YpontRadBarSx;
XpontRadDx = temp.XpontRadDx;
% YpontRadDx = temp.YpontRadDx;
% XpontRadSx = temp.XpontRadSx;
YpontRadSx = temp.YpontRadSx;
XpontSplitDx = temp.XpontSplitDx;
% YpontSplitDx = temp.YpontSplitDx;
% XpontSplitSx = temp.XpontSplitSx;
YpontSplitSx = temp.YpontSplitSx;
YpontSplitBarSx = temp.YpontSplitBarSx;

mOrto = temp.mOrto;

xc=[];
yc=[];      % (xc,yc)=PMs centers
xair=[];
yair=[];    % (xair,yair)=rotor air zone centers
xmag=[];
ymag=[];    % (xmag,ymag)=PMs magnetization direction
Br = [];
xmagair=[];
ymagair=[]; % (xmagair,ymagair)= air zone magnetization direction (maybe useful for FEMM definition???)

% fluxPortion=3;
for kk=1:nlay
    % PMs of the central segment
    xmedBar1=(B1k(kk)+B2k(kk))/2;
    ymedBar1=(temp.yvert2pt2(kk)+YpontRadBarDx(kk))/2;
    xc=[xc,xmedBar1];
    yc=[yc,ymedBar1];
    Br = [Br mat.LayerMag.Br(kk)];    % 1 block
    xmag=[xmag,cos(0)];
    ymag=[ymag,sin(0)];
    
    % Air zone at the end of the barrier (only if material~=Air and
    % Splitted internal ribs)
    if ((~strcmp(mat.LayerMag.MatName,'Air'))||(~isnan(XpontSplitDx(1,kk))))
        xtmp=(xxD1k(kk)+xxD2k(kk)+xpont(kk))/3;
        ytmp=(yyD1k(kk)+yyD2k(kk)+ypont(kk))/3;
        xair=[xair,xtmp];
        yair=[yair,ytmp];
    end
    %Br=[Br mat.LayerMag.Br(kk)]; % I'm not so sure...
    
    % Internal Radial Rib (only if material~=Air)
    if ((~strcmp(mat.LayerMag.MatName,'Air'))&&(~isnan(XpontRadDx(kk))))
        xtmp=(B1k(kk)+B2k(kk))/2;
        ytmp=(YpontRadBarSx(kk)+YpontRadSx(kk))/2;
        xair=[xair,xtmp];
        yair=[yair,ytmp];
    end
    
    % Splitted Internal Radial Rib (only if material~=Air)
    if (~strcmp(mat.LayerMag.MatName,'Air')&&~isnan(XpontSplitDx(1,kk)))
        xtmp=(B1k(kk)+B2k(kk))/2;
        ytmp=(YpontSplitBarSx(2,kk)+YpontSplitSx(2,kk))/2;
        xair=[xair,xtmp];
        yair=[yair,ytmp];
    end
    
    % Flux barrier arms (only if kk>1 && material~=Air)
    if ((kk>1)&&(~strcmp(mat.LayerMag.MatName,'Air'))&&(geo.Areavert(kk)>0))
        % PMs
        d2221 = calc_distanza_punti([XpBar2(kk),YpBar2(kk)],[xxD2k(kk),yyD2k(kk)]);
        if (d2221>0.5)
            [a_22tmp,b_22tmp,c_22tmp]=retta_per_2pti(temp.xob2pt1(kk),temp.yob2pt1(kk),temp.xob1pt2(kk),temp.yob1pt2(kk));
            [a_21tmp,b_21tmp,c_21tmp]=retta_per_2pti(temp.xob1pt1(kk),temp.yob1pt1(kk),temp.xob2pt2(kk),temp.yob2pt2(kk));
        else
            [a_22tmp,b_22tmp,c_22tmp]=retta_per_2pti(temp.xob2pt1(kk),temp.yob2pt1(kk),temp.xob1pt2(kk),temp.yob1pt2(kk));
            [a_21tmp,b_21tmp,c_21tmp]=retta_per_2pti(temp.xob1pt1(kk),temp.yob1pt1(kk),xpont(kk),ypont(kk));
            
        end
        [xmedBar3,ymedBar3]=intersezione_tra_rette(a_22tmp,b_22tmp,c_22tmp,a_21tmp,b_21tmp,c_21tmp);
        xc=[xc,xmedBar3];
        yc=[yc,ymedBar3];
        xmag=[xmag,cos(atan(mOrto))];
        ymag=[ymag,sin(atan(mOrto))];
        Br = [Br mat.LayerMag.Br(kk)];    % add another blocks for layers > 1
        
        % Air zone in the corner of the flux barrier
        if (YpBar2(kk)==yyD2k(kk))
            [a3,b3,c3]=retta_per_2pti(XpBar2(kk),YpBar2(kk),xpont(kk),ypont(kk));
        else
            [a3,b3,c3]=retta_per_2pti(XpBar2(kk),YpBar2(kk),xxD2k(kk),yyD2k(kk));
        end
        mOrto= b3/a3;
        
        xmedBar2=(XpBar1(kk)+XpBar2(kk))/2;
        ymedBar2=(YpMag1B1(kk)+YpBar1(kk))/2;
        [a2,b2,c2]=retta_per_2pti(XpMag1B1(kk),YpMag1B1(kk),XpBar2(kk),YpBar2(kk));
        mParal=-(a2/b2);
        %if strcmp (mat.LayerMag.MatName,'Bonded-Magnet') || strcmp(mat.LayerMag.MatName,'Air')
        xair=[xair,xmedBar2];
        yair=[yair,ymedBar2];
        xmagair=[xmagair,cos(atan(mParal))];
        ymagair=[ymagair,sin(atan(mParal))];
        Br = [Br mat.LayerMag.Br(kk)];
    end
end
zmag=zeros(1,size(xmag,2));

temp.xc=xc;
temp.yc=yc;
temp.xair=xair;
temp.yair=yair;
temp.xmag=xmag;
temp.ymag=ymag;
temp.xmagair=xmagair;
temp.ymagair=ymagair;
temp.zmag=zmag;

mat.LayerMag.Br = [Br Br];   % doubles Br pieces (half pole + half pole)





