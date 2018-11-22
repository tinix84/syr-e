
% Copyright 2018
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



% central point and direction of magnetization

function [geo,temp] = MagnetFullFill_Seg_fun(geo,mat,temp)

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
xpont = temp.xpont;
ypont = temp.ypont;
B1k = temp.B1k;
B2k = temp.B2k;

DTrasl = geo.DTrasl;

XpontRadBarDx = temp.XpontRadBarDx;
YpontRadBarDx = temp.YpontRadBarDx;
XpontRadBarSx = temp.XpontRadBarSx;
YpontRadBarSx = temp.YpontRadBarSx;
% XpontRadDx = temp.XpontRadDx;
% YpontRadDx = temp.YpontRadDx;
% XpontRadSx = temp.XpontRadSx;
% YpontRadSx = temp.YpontRadSx;
% XpontSplitDx = temp.XpontSplitDx;
% YpontSplitDx = temp.YpontSplitDx;
% XpontSplitSx = temp.XpontSplitSx;
% YpontSplitSx = temp.YpontSplitSx;
% YpontSplitBarSx = temp.YpontSplitBarSx;

% additional point
for kk=1:nlay
    if (yyD2k(kk)==YpBar2(kk))
        [a,b,c]=retta_per_2pti(XpBar2(kk),YpBar2(kk),xpont(kk),ypont(kk));
    else
        [a,b,c]=retta_per_2pti(XpBar2(kk),YpBar2(kk),xxD2k(kk),yyD2k(kk));
    end
    a2=b;
    b2=-a;
    c2=a*YpBar2(kk)-b*XpBar2(kk);
    [a3,b3,c3]=retta_per_2pti(XpBar1(kk),YpBar1(kk),xxD1k(kk),yyD1k(kk));
    [x,y]=intersezione_tra_rette(a2,b2,c2,a3,b3,c3);
    if x<XpBar1(kk)
        
        XpMag2B1(kk)=XpBar1(kk);
        YpMag2B1(kk)=YpBar1(kk);
        
    else
        XpMag2B1(kk)=x;
        YpMag2B1(kk)=y;
    end
    % internal barrier side intersection
    [a4,b4,c4]=retta_per_2pti(B2k(kk),0,XpBar2(kk),YpBar2(kk));
    [a5,b5,c5]=retta_per_2pti(B1k(kk),0,XpBar1(kk),YpBar1(kk));
    m6=b4/a4; d6=YpBar2(kk)-m6*XpBar2(kk);
    a6=-m6;
    b6=1;
    c6=-d6;
    [x,y]=intersezione_tra_rette(a5,b5,c5,a6,b6,c6);
    if (abs((YpBar1(kk)-y)+1j*(XpBar1(kk)-x))<0.5)
        XpMag1B1(kk)=XpBar1(kk);
        YpMag1B1(kk)=YpBar1(kk);
    elseif (y>YpBar1(kk))
        XpMag1B1(kk)=XpBar1(kk);
        YpMag1B1(kk)=YpBar1(kk);
    else
        XpMag1B1(kk)=x;
        YpMag1B1(kk)=y;
    end
end
% clear a b c a1 b1 c1 a2 b2 c2 a3 b3 c3 a4 b4 c4 x y

temp.XpMag2B1=XpMag2B1;
temp.YpMag2B1=YpMag2B1;
temp.XpMag1B1=XpMag1B1;
temp.YpMag1B1=YpMag1B1;

% evaluation of the direction of magnetization
xc=[];
yc=[];
xmag=[];
ymag=[]; Br = [];
fluxPortion=3;
for kk=1:nlay
    % #1 magnet portion per barrier
    [a1b,b1b,c1b]=retta_per_2pti(XpontRadBarSx(kk),YpontRadBarSx(kk),XpBar2(kk),YpBar2(kk));
    [a2b,b2b,c2b]=retta_per_2pti(XpontRadBarDx(kk),YpontRadBarDx(kk),XpMag1B1(kk),YpMag1B1(kk));
    
    %[xmedBar1,ymedBar1]=intersezione_tra_rette(a1b,b1b,c1b,a2b,b2b,c2b);
    xmedBar1=mean([XpontRadBarSx(kk),XpontRadBarDx(kk),XpBar1(kk),XpBar2(kk)]);
    ymedBar1=mean([YpontRadBarSx(kk),YpontRadBarDx(kk),YpBar1(kk),YpBar2(kk)]);
    xc=[xc,xmedBar1];
    yc=[yc,ymedBar1];
    Br = [Br mat.LayerMag.Br(kk)];    % 1 block
    
    [a1,b1,c1]=retta_per_2pti(B2k(kk),0,XpBar2(kk),YpBar2(kk));
    if (DTrasl==0)
        mOrto=0;
    else
        mOrto=b1/a1;
    end
    xmag=[xmag,cos(atan(mOrto))];
    ymag=[ymag,sin(atan(mOrto))];

    if (YpBar2(kk)==yyD2k(kk))
        [a3,b3,c3]=retta_per_2pti(XpBar2(kk),YpBar2(kk),xpont(kk),ypont(kk));
    else
        [a3,b3,c3]=retta_per_2pti(XpBar2(kk),YpBar2(kk),xxD2k(kk),yyD2k(kk));
    end
    mOrto=b3/a3;
    %xmedBar3=(XpMag2B1(kk)+xxD2k(kk))/2;
    %ymedBar3=(yyD1k(kk)+YpBar2(kk))/2;
    xmedBar3=mean([xxD1k(kk),xxD2k(kk),XpBar1(kk),XpBar2(kk)]);
    ymedBar3=mean([yyD1k(kk),yyD2k(kk),YpBar1(kk),YpBar2(kk)]);
    xc=[xc,xmedBar3];
    yc=[yc,ymedBar3];
    xmag=[xmag,cos(atan(mOrto))];
    ymag=[ymag,sin(atan(mOrto))];
    Br = [Br mat.LayerMag.Br(kk)];    % add another blocks to all layers
    
end

% zmag=zeros(1,size(xmag,2));

temp.mOrto = mOrto;



