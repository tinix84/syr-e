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

%% %%%%%%
%% 2013/11/04 MG file that assigns barrier points  
%% for ISeg geo and Seg geo
%% %%%%%%
%% i successivi export
%% Assegna_punti_banane_FluidBar
Bx0=temp.Bx0;
B1k=temp.B1k;
B2k=temp.B2k;
xpont=temp.xpont;
ypont=temp.ypont;
xxD1k=temp.xxD1k;
yyD1k=temp.yyD1k;
xxD2k=temp.xxD2k;
yyD2k=temp.yyD2k;
XpBar1=temp.XpBar1;
YpBar1=temp.YpBar1;
XpBar2=temp.XpBar2;
YpBar2=temp.YpBar2;

xTraf1=temp.xTraf1;
xTraf2=temp.xTraf2;
yTraf1=temp.yTraf1;
yTraf2=temp.yTraf2;

arcLayTraf1=temp.arcLayTraf1;
arcLayTraf2=temp.arcLayTraf2;

XcRibTraf1=temp.XcRibTraf1;
YcRibTraf1=temp.YcRibTraf1;
XcRibTraf2=temp.XcRibTraf2;
YcRibTraf2=temp.YcRibTraf2;

XpontRadSx=temp.XpontRadSx;
YpontRadSx=temp.YpontRadSx;
XpontRadDx=temp.XpontRadDx;
YpontRadDx=temp.YpontRadDx;
XpontRadBarDx=temp.XpontRadBarDx;
YpontRadBarDx=temp.YpontRadBarDx;
XpontRadBarSx=temp.XpontRadBarSx;
YpontRadBarSx=temp.YpontRadBarSx;
%% Next 4 statements are redundant but ensures compatibility with script
%% that assigns points at circular barriers
XBan1dx=temp.B2k;  
YBan1dx=temp.YpontRadBarDx;
XBan1sx=temp.B1k;    
YBan1sx=temp.YpontRadBarSx;
%% Error mex no linear barrier boundary
error_mex=temp.error_mex;

xc=temp.xc;
yc=temp.yc;

rotore=[];

for ii=1:geo.nlay
    
    if (YpontRadSx(ii)~=0)
        rotore=[rotore;XpontRadSx(ii),YpontRadSx(ii),XpontRadDx(ii),YpontRadDx(ii),NaN,NaN,0;...
            XpontRadSx(ii),YpontRadSx(ii),XpontRadBarSx(ii) YpontRadBarSx(ii),NaN,NaN,0;...
            XpontRadDx(ii),YpontRadDx(ii),XpontRadBarDx(ii) YpontRadBarDx(ii),NaN,NaN,0];
        
    end
    if (error_mex(ii)==0)
        rotore=[rotore;
            XpontRadBarSx(ii) YpontRadBarSx(ii) XpBar1(ii) YpBar1(ii) NaN NaN 0;
            XpontRadBarDx(ii) YpontRadBarDx(ii) XpBar2(ii) YpBar2(ii) NaN NaN 0;
            XpBar1(ii) YpBar1(ii) xxD1k(ii) yyD1k(ii) NaN NaN 0;
            XpBar2(ii) YpBar2(ii) xxD2k(ii) yyD2k(ii) NaN NaN 0];
        rotore=[ rotore;
            XcRibTraf1(ii) YcRibTraf1(ii) xpont(ii) ypont(ii) xxD1k(ii) yyD1k(ii) 1];
        rotore=[rotore;
            XcRibTraf2(ii) YcRibTraf2(ii) xxD2k(ii) yyD2k(ii) xpont(ii) ypont(ii) 1];
    else
        rotore=[ rotore;
            XcRibTraf1(ii) YcRibTraf1(ii) xxD2k(ii) yyD2k(ii) xxD1k(ii) yyD1k(ii) 1];
        rotore=[rotore;
            XcRibTraf2(ii) YcRibTraf2(ii) xxD1k(ii) yyD1k(ii) xxD2k(ii) yyD2k(ii) 1];
        
    end
    
end

clear temp
