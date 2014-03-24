% Determinazione del centro e raggio di una circonferenza per 3 pti noti

function [xc,yc,r,tA,tB]=circonferenza_per_3_pti(xA,yA,xB,yB,xC,yC)

delta=det([xA,yA,1; xB,yB,1;xC,yC,1]);
deltaA=det([-(xA^2+yA^2),yA,1;-(xB^2+yB^2),yB,1;-(xC^2+yC^2),yC,1]);
deltaB=det([xA,-(xA^2+yA^2),1;xB,-(xB^2+yB^2),1;xC,-(xC^2+yC^2),1]);
deltaC=det([xA,yA,-(xA^2+yA^2);xB,yB,-(xB^2+yB^2);xC,yC,-(xC^2+yC^2)]);
a=deltaA/delta; b=deltaB/delta; c=deltaC/delta;
xc=-a/2;
yc=-b/2;
r=sqrt(xc^2+yc^2-c);
tA=atan2((yA-yc),(xA-xc));
tB=atan2((yB-yc),(xB-xc));
% angle=(tA-tB)*180/pi;

end