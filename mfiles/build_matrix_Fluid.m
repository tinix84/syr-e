%% %%%%%%
%% 2013/07/05 MG file di assegnazione dei punti delle barriere di flux per
%% %%%%%%
%% i successivi export
%% Assegna_punti_banane_FluidBar

% Riassegnazione punti centrali di barriera
xc=temp.xpont1;
yc=temp.ypont1;
XBan1sx=temp.B1k;
%%%%%%%%%%%%%%%%%%%

Bx0=temp.Bx0;
By0=temp.By0;
B1k=temp.B1k;
B2k=temp.B2k;
xpont=temp.xpont;
ypont=temp.ypont;

xpont1=temp.xpont1;
ypont1=temp.ypont1;

xTraf1=temp.xTraf1;
xTraf2=temp.xTraf2;
yTraf1=temp.yTraf1;
yTraf2=temp.yTraf2;
arcLayer1=temp.arcLayer1;
arcLayer2=temp.arcLayer2;

arcLayTraf1=temp.arcLayTraf1;
arcLayTraf2=temp.arcLayTraf2;

XcRibTraf1=temp.XcRibTraf1;
YcRibTraf1=temp.YcRibTraf1;
XcRibTraf2=temp.XcRibTraf2;
YcRibTraf2=temp.YcRibTraf2;
%% Punti per i ribs radiali
XpontRadSx=temp.XpontRadSx;
YpontRadSx=temp.YpontRadSx;
XpontRadDx=temp.XpontRadDx;
YpontRadDx=temp.YpontRadDx;
XpontRadBarDx=temp.XpontRadBarDx;
YpontRadBarDx=temp.YpontRadBarDx;
XpontRadBarSx=temp.XpontRadBarSx;
YpontRadBarSx=temp.YpontRadBarSx;

xxD1k=temp.xxD1k;
yyD1k=temp.yyD1k;
xxD2k=temp.xxD2k;
yyD2k=temp.yyD2k;

XcBar1=temp.XcBar1;
YcBar1=temp.YcBar1;
XcBar2=temp.XcBar2;
YcBar2=temp.YcBar2;
XcBarLast_mean=temp.XcBarLast_mean;
YcBarLast_mean=temp.YcBarLast_mean;

xxB1k_mean=temp.xxB1k_mean;
yyB1k_mean=temp.yyB1k_mean;
xxB2k_mean=temp.xxB2k_mean;
yyB2k_mean=temp.yyB2k_mean;
xxB1k_mean3=temp.xxB1k_mean3;
yyB1k_mean3=temp.yyB1k_mean3;
xxB1k_mean2=temp.xxB1k_mean2;
yyB1k_mean2=temp.yyB1k_mean2;

arcLayer3=temp.arcLayer3;

rotore=[];

for ii=1:geo.nlay
    
    if (YpontRadSx(ii)~=0)
        rotore=[rotore;XpontRadSx(ii),YpontRadSx(ii),XpontRadDx(ii),YpontRadDx(ii),NaN,NaN,0;...
            XpontRadSx(ii),YpontRadSx(ii),XpontRadBarSx(ii) YpontRadBarSx(ii),NaN,NaN,0;...
            XpontRadDx(ii),YpontRadDx(ii),XpontRadBarDx(ii) YpontRadBarDx(ii),NaN,NaN,0];
        
        rotore=[ rotore;
            XcRibTraf1(ii) YcRibTraf1(ii) xpont(ii) ypont(ii) xxD1k(ii) yyD1k(ii) 1];
        rotore=[rotore;
            XcRibTraf2(ii) YcRibTraf2(ii) xxD2k(ii) yyD2k(ii) xpont(ii) ypont(ii) 1];
        rotore=[rotore;
            XcBar2(ii) YcBar2(ii) xxD2k(ii) yyD2k(ii) XpontRadBarDx(ii) YpontRadBarDx(ii) 1];
        if (ii<nlay || nlay==1)
            rotore=[rotore;
                XcBar1(ii) YcBar1(ii) xxD1k(ii) yyD1k(ii) XpontRadBarSx(ii) YpontRadBarSx(ii) 1];
        else
            rotore=[rotore;
                XpontRadBarSx(ii) YpontRadBarSx(ii) xxB1k_mean2 yyB1k_mean2 NaN NaN 0;
                xxB1k_mean3 yyB1k_mean3 xxD1k(ii) yyD1k(ii) NaN NaN 0;
                XcBarLast_mean YcBarLast_mean xxB1k_mean3 yyB1k_mean3 xxB1k_mean2 yyB1k_mean2 1];
        end
        
    else
        rotore=[ rotore;
            XcRibTraf1(ii) YcRibTraf1(ii) xpont(ii) ypont(ii) xxD1k(ii) yyD1k(ii) 1];
        rotore=[rotore;
            XcRibTraf2(ii) YcRibTraf2(ii) xxD2k(ii) yyD2k(ii) xpont(ii) ypont(ii) 1];
        rotore=[rotore;
            XcBar2(ii) YcBar2(ii) xxD2k(ii) yyD2k(ii) XpontRadBarDx(ii) YpontRadBarDx(ii) 1];
        if (ii<nlay || nlay==1)
            rotore=[rotore;
                XcBar1(ii) YcBar1(ii) xxD1k(ii) yyD1k(ii) XpontRadBarSx(ii) YpontRadBarSx(ii) 1];
        else
            rotore=[rotore;
                XpontRadBarSx(ii) YpontRadBarSx(ii) xxB1k_mean2 yyB1k_mean2 NaN NaN 0;
                xxB1k_mean3 yyB1k_mean3 xxD1k(ii) yyD1k(ii) NaN NaN 0;
                XcBarLast_mean YcBarLast_mean xxB1k_mean3 yyB1k_mean3 xxB1k_mean2 yyB1k_mean2 1];
        end
    end
end
