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

function rotore = build_matrix_Seg(temp,geo)

xc=temp.xc;
yc=temp.yc;
XBan1sx=temp.B1k;

Bx0=temp.Bx0;
B1k=temp.B1k;
B2k=temp.B2k;
xpont=temp.xpont;
ypont=temp.ypont;

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


XpBar1=temp.XpBar1;
YpBar1=temp.YpBar1;
XpBar2=temp.XpBar2;
YpBar2=temp.YpBar2;

xxD1k=temp.xxD1k;
yyD1k=temp.yyD1k;
xxD2k=temp.xxD2k;
yyD2k=temp.yyD2k;

% Coordinate per raccordi barriere di flux:
XcRacc_B1=temp.XcRacc_B1;
YcRacc_B1=temp.YcRacc_B1;
xRaccR1_B1=temp.xRaccR1_B1;
yRaccR1_B1=temp.yRaccR1_B1;
xRaccR2_B1=temp.xRaccR2_B1;
yRaccR2_B1=temp.yRaccR2_B1;
%
XcRacc_B2=temp.XcRacc_B2;
YcRacc_B2=temp.YcRacc_B2;
xRaccR1_B2=temp.xRaccR1_B2;
yRaccR1_B2=temp.yRaccR1_B2;
xRaccR2_B2=temp.xRaccR2_B2;
yRaccR2_B2=temp.yRaccR2_B2;

R_RaccB1=temp.R_RaccB1;
R_RaccB2=temp.R_RaccB2;

%% Punti per i ribs radiali
XpontRadSx=temp.XpontRadSx;
YpontRadSx=temp.YpontRadSx;
XpontRadDx=temp.XpontRadDx;
YpontRadDx=temp.YpontRadDx;
XpontRadBarDx=temp.XpontRadBarDx;
YpontRadBarDx=temp.YpontRadBarDx;
XpontRadBarSx=temp.XpontRadBarSx;
YpontRadBarSx=temp.YpontRadBarSx;
%% Additional point for magnet
% if (geo.Br~=0)
XpMag2B1=temp.XpMag2B1;
YpMag2B1=temp.YpMag2B1;
XpMag1B1=temp.XpMag1B1;
YpMag1B1=temp.YpMag1B1;
% end
%%
LowDimBarrier=temp.LowDimBarrier;
%%
rotore=[];


for ii=1:geo.nlay
    
    if (YpontRadSx(ii)~=0)
%         if (ii==1 && LowDimBarrier(ii)==1)
            %             rotore=[rotore;(B1k(ii)+B2k(ii))/2,0,B1k(ii),0,B2k(ii),0,1;
            %                 (B1k(ii)+B2k(ii))/2,0,B2k(ii),0,B1k(ii),0,1];
%         else
            rotore=[rotore;XpontRadSx(ii),YpontRadSx(ii),XpontRadDx(ii),YpontRadDx(ii),NaN,NaN,0;...
                XpontRadSx(ii),YpontRadSx(ii),XpontRadBarSx(ii) YpontRadBarSx(ii),NaN,NaN,0;...
                XpontRadDx(ii),YpontRadDx(ii),XpontRadBarDx(ii) YpontRadBarDx(ii),NaN,NaN,0];
            
            if (R_RaccB1(ii)<=0.5)
                rotore=[rotore;
                    XpontRadBarSx(ii) YpontRadBarSx(ii) XpBar1(ii) YpBar1(ii) NaN NaN 0;
                    XpBar1(ii) YpBar1(ii) xxD1k(ii) yyD1k(ii) NaN NaN 0];
            else
                rotore=[rotore;
                    XcRacc_B1(ii) YcRacc_B1(ii) xRaccR2_B1(ii) yRaccR2_B1(ii) xRaccR1_B1(ii) yRaccR1_B1(ii) 1;
                    XpontRadBarSx(ii) YpontRadBarSx(ii) xRaccR1_B1(ii) yRaccR1_B1(ii) NaN NaN 0;
                    xRaccR2_B1(ii) yRaccR2_B1(ii) xxD1k(ii) yyD1k(ii) NaN NaN 0];
            end
            
            if (R_RaccB2(ii)<=0.5)
                rotore=[rotore;
                    XpontRadBarDx(ii) YpontRadBarDx(ii) XpBar2(ii) YpBar2(ii) NaN NaN 0;
                    XpBar2(ii) YpBar2(ii) xxD2k(ii) yyD2k(ii) NaN NaN 0];
            else
                rotore=[rotore;
                    XcRacc_B2(ii) YcRacc_B2(ii) xRaccR2_B2(ii) yRaccR2_B2(ii) xRaccR1_B2(ii) yRaccR1_B2(ii) 1;
                    XpontRadBarDx(ii) YpontRadBarDx(ii) xRaccR1_B2(ii) yRaccR1_B2(ii) NaN NaN 0;
                    xRaccR2_B2(ii) yRaccR2_B2(ii) xxD2k(ii) yyD2k(ii) NaN NaN 0];
                
            end
            rotore=[ rotore;
                XcRibTraf1(ii) YcRibTraf1(ii) xpont(ii) ypont(ii) xxD1k(ii) yyD1k(ii) 1];
            rotore=[rotore;
                XcRibTraf2(ii) YcRibTraf2(ii) xxD2k(ii) yyD2k(ii) xpont(ii) ypont(ii) 1];
            
%         end
    else
%         if (ii==1 && LowDimBarrier(ii)==1)
            %             rotore=[rotore;B1k(ii),0,(B1k(ii)+B2k(ii))/2,ypont(1),NaN,NaN,0;
            %                 (B1k(ii)+B2k(ii))/2,ypont(1),B2k(ii),0,NaN,NaN,0];
%         else
            
            if (R_RaccB1(ii)<=0.5)
                rotore=[rotore;
                    B1k(ii) 0 XpBar1(ii) YpBar1(ii) NaN NaN 0;
                    XpBar1(ii) YpBar1(ii) xxD1k(ii) yyD1k(ii) NaN NaN 0];
            else
                rotore=[rotore;
                    XcRacc_B1(ii) YcRacc_B1(ii) xRaccR2_B1(ii) yRaccR2_B1(ii) xRaccR1_B1(ii) yRaccR1_B1(ii) 1;
                    B1k(ii) 0 xRaccR1_B1(ii) yRaccR1_B1(ii) NaN NaN 0;
                    xRaccR2_B1(ii) yRaccR2_B1(ii) xxD1k(ii) yyD1k(ii) NaN NaN 0];
            end
            
            if (R_RaccB2(ii)<=0.5)
                rotore=[rotore;
                    B2k(ii) 0 XpBar2(ii) YpBar2(ii) NaN NaN 0;
                    XpBar2(ii) YpBar2(ii) xxD2k(ii) yyD2k(ii) NaN NaN 0];
            else
                rotore=[rotore;
                    XcRacc_B2(ii) YcRacc_B2(ii) xRaccR2_B2(ii) yRaccR2_B2(ii) xRaccR1_B2(ii) yRaccR1_B2(ii) 1;
                    B2k(ii) 0 xRaccR1_B2(ii) yRaccR1_B2(ii) NaN NaN 0;
                    xRaccR2_B2(ii) yRaccR2_B2(ii) xxD2k(ii) yyD2k(ii) NaN NaN 0];
                
            end
            rotore=[ rotore;
                XcRibTraf1(ii) YcRibTraf1(ii) xpont(ii) ypont(ii) xxD1k(ii) yyD1k(ii) 1];
            rotore=[rotore;
                XcRibTraf2(ii) YcRibTraf2(ii) xxD2k(ii) yyD2k(ii) xpont(ii) ypont(ii) 1];
%         end
    end
    
    %     if (geo.Br(ii)~=0)
    if (YpontRadSx(ii)==0)
        rotore=[rotore;
            B1k(ii),0,B2k(ii),0,NaN,NaN,0];
    end
    
    if (R_RaccB2(ii)>0.5 &&R_RaccB1(ii)<=0.5)
        rotore=[rotore;
            xRaccR2_B2(ii) yRaccR2_B2(ii) XpMag2B1(ii) YpMag2B1(ii) NaN NaN 0;
            xRaccR1_B2(ii) yRaccR1_B2(ii) XpMag1B1(ii) YpMag1B1(ii) NaN NaN 0];
        
    elseif (R_RaccB1(ii)>0.5 && R_RaccB2(ii)<=0.5)
        rotore=[rotore;
            xRaccR2_B1(ii) yRaccR2_B1(ii) XpMag2B1(ii) YpMag2B1(ii) NaN NaN 0;
            xRaccR1_B1(ii) yRaccR1_B1(ii) XpMag1B1(ii) YpMag1B1(ii) NaN NaN 0];
        
    elseif (R_RaccB1(ii)>0.5 && R_RaccB2(ii)>0.5)
        rotore=[rotore;
            xRaccR2_B2(ii) yRaccR2_B2(ii) xRaccR2_B1(ii) yRaccR2_B1(ii) NaN NaN 0;
            xRaccR1_B2(ii) yRaccR1_B2(ii) xRaccR1_B1(ii) yRaccR1_B1(ii) NaN NaN 0];
        
    else
        %  segment for the barrier
%         if ii == 1
%             rotore=[rotore;
%                 XpBar1(ii)  YpBar1(ii) XpMag1B1(ii) YpMag1B1(ii) NaN NaN 0];
%         else
            rotore=[rotore;
                XpBar2(ii)  YpBar2(ii) XpMag1B1(ii) YpMag1B1(ii) NaN NaN 0];
%         end
    end
end


