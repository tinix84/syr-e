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

function rotore = build_matrix_ISeg(temp,geo)

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
% radial ribs coordinate
XpontRadSx=temp.XpontRadSx;
YpontRadSx=temp.YpontRadSx;
XpontRadDx=temp.XpontRadDx;
YpontRadDx=temp.YpontRadDx;
XpontRadBarDx=temp.XpontRadBarDx;
YpontRadBarDx=temp.YpontRadBarDx;
XpontRadBarSx=temp.XpontRadBarSx;
YpontRadBarSx=temp.YpontRadBarSx;
% Additional point for magnet
XpMag1B1=temp.XpMag1B1;
YpMag1B1=temp.YpMag1B1;

% Error mex no linear barrier boundary
error_mex=temp.error_mex;
% elementi per completamento punti
Mag=temp.Mag;

xob1pt1=temp.xob1pt1;
yob1pt1=temp.yob1pt1;
xob1pt2=temp.xob1pt2;
yob1pt2=temp.yob1pt2;
xob2pt1=temp.xob2pt1;
yob2pt1=temp.yob2pt1;
xob2pt2=temp.xob2pt2;
yob2pt2=temp.yob2pt2;

   
rotore=[];

for ii=1:geo.nlay
    
    if (YpontRadSx(ii)~=0)
        rotore=[rotore;XpontRadSx(ii),YpontRadSx(ii),XpontRadDx(ii),YpontRadDx(ii),NaN,NaN,0;...
            XpontRadSx(ii),YpontRadSx(ii),XpontRadBarSx(ii) YpontRadBarSx(ii),NaN,NaN,0;...
            XpontRadDx(ii),YpontRadDx(ii),XpontRadBarDx(ii) YpontRadBarDx(ii),NaN,NaN,0];
    else
        rotore=[rotore;
            B1k(ii),0,B2k(ii),0,NaN,NaN,0];
    end
    if (error_mex(ii)==0)
        if ii==1
            rotore=[rotore;
                XpontRadBarSx(ii) YpontRadBarSx(ii) XpBar1(ii) YpBar1(ii) NaN NaN 0;
                XpontRadBarDx(ii) YpontRadBarDx(ii) XpBar2(ii) YpBar2(ii) NaN NaN 0;
                XpBar1(ii) YpBar1(ii) xxD1k(ii) yyD1k(ii) NaN NaN 0;
                XpBar2(ii) YpBar2(ii) xxD2k(ii) yyD2k(ii) NaN NaN 0];
            rotore=[rotore;
                XcRibTraf1(ii) YcRibTraf1(ii) xpont(ii) ypont(ii) xxD1k(ii) yyD1k(ii) 1];
            rotore=[rotore;
                XcRibTraf2(ii) YcRibTraf2(ii) xxD2k(ii) yyD2k(ii) xpont(ii) ypont(ii) 1];
        else
                rotore=[rotore;
                    XpontRadBarSx(ii) YpontRadBarSx(ii) XpBar1(ii) YpBar1(ii) NaN NaN 0;
                    XpontRadBarDx(ii) YpontRadBarDx(ii) XpBar2(ii) YpBar2(ii) NaN NaN 0;
                    XpBar1(ii) YpBar1(ii) xob1pt1(ii) yob1pt1(ii) NaN NaN 0;
                    xob1pt1(ii) yob1pt1(ii) xob1pt2(ii) yob1pt2(ii) NaN NaN 0;
                    xob1pt2(ii) yob1pt2(ii) xxD1k(ii) yyD1k(ii) NaN NaN 0;
                    XpBar2(ii) YpBar2(ii) xob2pt2(ii) yob2pt2(ii) NaN NaN 0;
                    xob2pt2(ii) yob2pt2(ii) xxD2k(ii) yyD2k(ii) NaN NaN 0];
                rotore=[rotore;
                    XcRibTraf1(ii) YcRibTraf1(ii) xpont(ii) ypont(ii) xxD1k(ii) yyD1k(ii) 1];
                rotore=[rotore;
                    XcRibTraf2(ii) YcRibTraf2(ii) xxD2k(ii) yyD2k(ii) xpont(ii) ypont(ii) 1];

        end
        
    else
        rotore=[rotore;
            XcRibTraf1(ii) YcRibTraf1(ii) xxD2k(ii) yyD2k(ii) xxD1k(ii) yyD1k(ii) 1];
        rotore=[rotore;
            XcRibTraf2(ii) YcRibTraf2(ii) xxD1k(ii) yyD1k(ii) xxD2k(ii) yyD2k(ii) 1];
        
    end
    

    
end

%% Inserimento matrice Mag

rotore=[rotore;Mag];

end

