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

function rotore = build_matrix_Circ_dx(temp,geo)

x0=geo.x0;

% XBan1dx=temp.XBanqdx;
% XBan1sx=temp.XBanqsx;
% xc=temp.xc;
% yc=temp.yc;
% XBan3dx=temp.X3;
% YBan3dx=temp.Y3;
% XBan3sx=temp.X4;
% YBan3sx=temp.Y4;
% 
% XpontRadDx=temp.XpontRadDx;
% YpontRadDx=temp.YpontRadDx;
% XpontRadSx=temp.XpontRadSx;
% YpontRadSx=temp.YpontRadSx;
% XpontRadBarDx=temp.XpontRadBarDx;
% YpontRadBarDx=temp.YpontRadBarDx;
% XpontRadBarSx=temp.XpontRadBarSx;
% YpontRadBarSx=temp.YpontRadBarSx;
% 
% YBan1dx=temp.YpontRadBarDx;
% YBan1sx=temp.YpontRadBarSx;
% error_mex=temp.error_mex;

B1k=temp.B1k;
B2k=temp.B2k;
xpont=temp.xpont;
ypont=temp.ypont;
xxD1k=temp.xxD1k;
yyD1k=temp.yyD1k;
xxD2k=temp.xxD2k;
yyD2k=temp.yyD2k;
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

Xmag5dx=temp.X5;
Ymag5dx=temp.Y5;
Xmag6sx=temp.X6;
Ymag6sx=temp.Y6;

xHalfDx=temp.xHalfDx;
yHalfDx=temp.yHalfDx;
xHalfSx=temp.xHalfSx;
yHalfSx=temp.yHalfSx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (geo.BarFillFac==0)
    rotore=[];
    
    for ii=1:geo.nlay
        
        if (YpontRadSx(ii)~=0) % ponticello radiale
            rotore=[rotore;
                XpontRadSx(ii),YpontRadSx(ii),XpontRadDx(ii),YpontRadDx(ii),NaN,NaN,0;
                XpontRadSx(ii),YpontRadSx(ii),XpontRadBarSx(ii),YpontRadBarSx(ii),NaN,NaN,0;
                XpontRadDx(ii),YpontRadDx(ii),XpontRadBarDx(ii),YpontRadBarDx(ii),NaN,NaN,0];
        else
            rotore=[rotore;
                XpontRadBarSx(ii), YpontRadBarSx(ii), XpontRadBarDx(ii), YpontRadBarDx(ii), NaN, NaN, 0];
        end
        % archi al traferro
        rotore=[ rotore;
            XcRibTraf1(ii) YcRibTraf1(ii) xpont(ii) ypont(ii) xxD1k(ii) yyD1k(ii) 1;
            XcRibTraf2(ii) YcRibTraf2(ii) xxD2k(ii) yyD2k(ii) xpont(ii) ypont(ii) 1];
        % archi laterali barriere
        rotore=[rotore;
            x0, 0, xxD1k(ii), yyD1k(ii), XpontRadBarSx(ii), YpontRadBarSx(ii), 1;
            x0, 0, xxD2k(ii), yyD2k(ii), XpontRadBarDx(ii), YpontRadBarDx(ii), 1];
        % 
        
%         if(error_mex(ii)==0) % disegno archi laterali e finali barriera
%             rotore=[ rotore;
%                 xc(ii) yc(ii) XBan3dx(ii) YBan3dx(ii) XBan3sx(ii) YBan3sx(ii) 1;
%                 x0 0 XBan3sx(ii) YBan3sx(ii) XBan1sx(ii) YBan1sx(ii) 1];
%             rotore=[rotore;
%                 x0 0 XBan3dx(ii) YBan3dx(ii) XBan1dx(ii) YBan1dx(ii) 1];
%         else
%             rotore=[rotore;
%                 xc(ii) 0 XBan3dx(ii) YBan3dx(ii) XBan3sx(ii) YBan3sx(ii) 1];
%             
%         end
        
        
        if (YpontRadSx(ii)==0) % bisettrice polo in barriera
            rotore=[rotore;
                xxD1k(ii) yyD1k(ii) xxD1k(ii) yyD1k(ii) NaN NaN 0];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    
    rotore=[];
    
    for ii=1:geo.nlay
        
        if (YpontRadSx(ii)~=0) % ponticello radiale
            rotore=[rotore;
                XpontRadSx(ii),YpontRadSx(ii),XpontRadDx(ii),YpontRadDx(ii),NaN,NaN,0;
                XpontRadSx(ii),YpontRadSx(ii),XpontRadBarSx(ii),YpontRadBarSx(ii),NaN,NaN,0;
                XpontRadDx(ii),YpontRadDx(ii),XpontRadBarDx(ii),YpontRadBarDx(ii),NaN,NaN,0];
        else
            rotore=[rotore;
                XpontRadBarSx(ii), YpontRadBarSx(ii), XpontRadBarDx(ii), YpontRadBarDx(ii), NaN, NaN, 0];
        end
        % archi al traferro
        rotore=[rotore;
            XcRibTraf1(ii) YcRibTraf1(ii) xpont(ii) ypont(ii) xxD1k(ii) yyD1k(ii) 1;
            XcRibTraf2(ii) YcRibTraf2(ii) xxD2k(ii) yyD2k(ii) xpont(ii) ypont(ii) 1];
        
        % archi barriere
        % rotore=[rotore;
        %     x0, 0, xxD1k(ii), yyD1k(ii), Xmag6sx(ii), Ymag6sx(ii), 1;
        %     x0, 0, Xmag6sx(ii), Ymag6sx(ii), XpontRadBarSx(ii), YpontRadBarSx(ii), 1;
        %     x0, 0, xxD2k(ii), yyD2k(ii), Xmag5dx(ii), Ymag5dx(ii), 1;
        %     x0, 0, Xmag5dx(ii), Ymag5dx(ii), XpontRadBarDx(ii), YpontRadBarDx(ii), 1
        %     Xmag6sx(ii) Ymag6sx(ii) Xmag5dx(ii) Ymag5dx(ii) NaN NaN 0;];
        
        rotore=[rotore;
            x0, 0, xxD1k(ii), yyD1k(ii), XpontRadBarSx(ii), YpontRadBarSx(ii), 1;
            x0, 0, xxD2k(ii), yyD2k(ii), XpontRadBarDx(ii), YpontRadBarSx(ii), 1];
                
        % bar fill factor e segmenti magnete
        rotore=[rotore;
            Xmag6sx(ii), Ymag6sx(ii), Xmag5dx(ii), Ymag5dx(ii), NaN, NaN, 0;
            xHalfSx(ii), yHalfSx(ii), xHalfDx(ii), yHalfDx(ii), NaN, NaN, 0];
        
%         %% Connect arcs one by one to avoid disconnection conflict
%         rotore=[rotore;
%             x0, 0, xxD2k(ii), yyD2k(ii),Xmag5dx(ii), Ymag5dx(ii),  1;
%             x0, 0, xxD1k(ii), yyD1k(ii),Xmag6sx(ii), Ymag6sx(ii), 1;
%             x0, 0, Xmag5dx(ii), Ymag5dx(ii), xHalfDx(ii), yHalfDx(ii), 1;
%             x0, 0, Xmag6sx(ii), Ymag6sx(ii), xHalfSx(ii), yHalfSx(ii), 1;
%             x0, 0, xHalfDx(ii), yHalfDx(ii), XpontRadBarDx(ii),YpontRadBarDx(ii), 1;
%             x0, 0, xHalfSx(ii), yHalfSx(ii), XpontRadBarSx(ii),YpontRadBarSx(ii), 1;];
%         
%         rotore=[rotore;
%             Xmag6sx(ii), Ymag6sx(ii), Xmag5dx(ii), Ymag5dx(ii), NaN, NaN, 0;
%             xHalfSx(ii), yHalfSx(ii), xHalfDx(ii), yHalfDx(ii), NaN, NaN, 0];
        %% 
%         if(error_mex(ii)==0)
%             rotore =[ rotore;
%                 xc(ii) yc(ii) XBan3dx(ii) YBan3dx(ii) XBan3sx(ii) YBan3sx(ii) 1;
%                 x0 0 XBan3sx(ii) YBan3sx(ii) Xmag6sx(ii) Ymag6sx(ii) 1;
%                 x0 0 Xmag6sx(ii) Ymag6sx(ii)  XBan1sx(ii) YBan1sx(ii) 1];
%             
%             rotore=[rotore;
%                 x0 0 XBan3dx(ii) YBan3dx(ii) Xmag5dx(ii) Ymag5dx(ii) 1;
%                 x0 0 Xmag5dx(ii) Ymag5dx(ii) XBan1dx(ii) YBan1dx(ii) 1];
%             
%         else
%             rotore=[rotore;
%                 xc(ii) 0 XBan3dx(ii) YBan3dx(ii) XBan3sx(ii) YBan3sx(ii) 1];
%             
%         end
        
%         if (YpontRadSx(ii)==0)
%             rotore=[rotore;
%                 XpontRadSx(ii) YpontRadSx(ii) XpontRadDx(ii) YpontRadDx(ii) NaN NaN 0];
%         end
    end
end

