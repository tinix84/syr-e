%% 2013/07/26 MG
%% Costruisce la matrice di linee ed archi al traferro, sono comprese nella construzione le linee ed archi di statore e di rotore...
% [traferro]: uscita viene data una matrice di archi e linee per la costruzione del
% traferro
% CENTRItraf.xy: matrice di centri dei materiali
function [traferro,CENTRItraf]=traferroMatr(pc,xr,g,Q,Qs,ps,p,meshTraf)

traferro=[];

%% Punti per il disegno delle linee ed archi al traferro:

[xArcStat1,yArcStat1] = rot_point(xr+2/3*g,0,-pc*pi/180);
[xArcStat2,yArcStat2] = rot_point(xr+2/3*g,0,(2*Qs-1)*pc*pi/180);

[xs1,ys1]=rot_point(xr+g,0,-pc*pi/180);
[xs2,ys2]=rot_point(xr+g,0,(2*Qs-1)*pc*pi/180);

[xArcTraf1,yArcTraf1] = rot_point(xr+1/3*g,0,-pc*pi/180);
[xArcTraf2,yArcTraf2] = rot_point(xr+1/3*g,0,(2*Qs-1)*pc*pi/180);

[xArcRot1,yArcRot1] = rot_point(xr+1/3*g,0,0);
[xArcRot2,yArcRot2] = rot_point(xr+1/3*g,0,ps*180/p*pi/180);
[xR1,yR1] = rot_point(xr,0,0);
[xR2,yR2] = rot_point(xr,0,ps*180/p*pi/180);

traferro=[xs1,ys1,xArcStat1,yArcStat1,NaN,NaN,0; 
          xArcStat1,yArcStat1,xArcTraf1,yArcTraf1,NaN,NaN,0;
          xs2,ys2,xArcStat2,yArcStat2,NaN,NaN,0; 
          xArcStat2,yArcStat2,xArcTraf2,yArcTraf2,NaN,NaN,0;
          xR1,yR1,xArcRot1,yArcRot1,NaN,NaN,0; 
          xR2,yR2,xArcRot2,yArcRot2,NaN,NaN,0; 
          0 0 xArcStat1 yArcStat1 xArcStat2 yArcStat2 1;
          0 0 xArcTraf1 yArcTraf1 xArcTraf2 yArcTraf2 1;
          0 0 xArcRot1 yArcRot1 xArcRot2 yArcRot2 1;
          ];

%% Punti Assegnazione dei materiali:
[xAirTrafSt,yAirTrafSt] = rot_point(xr+5/6*g,0,0.5*ps*180/p*pi/180);
[xAirTrafTr,yAirTrafTr] = rot_point(xr+3/6*g,0,0.5*ps*180/p*pi/180);
[xAirTrafRot,yAirTrafRot] = rot_point(xr+1/6*g,0,0.5*ps*180/p*pi/180);
CENTRItraf.xy=[xAirTrafSt,yAirTrafSt,2,meshTraf,0;
                xAirTrafTr,yAirTrafTr,2,meshTraf,0;
                xAirTrafRot,yAirTrafRot,2,meshTraf,1];

CENTRItraf.names={'air_gap_stator1','air_gap_stator2','air_gap_rotor'}';

end
