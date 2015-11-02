% disegna la linea di confine statore - rotore che va rifatta ad ogni
% posizione

function disegna_bordo_mobile(geo,th_m)

xr = geo.xr;
g = geo.g;
p = geo.p;
ns = geo.ns;

pc = 360/(ns*p)/2;  % mezzo passo cava

group = 20;

x0 = xr+1/3*g; y0 = 0;

[x1,y1] = rot_point(x0,y0,-pc*pi/180);
[x2,y2] = rot_point(x0,y0,th_m*pi/180);
[x3,y3] = rot_point(x0,y0,(-pc + 180/p)*pi/180);
[x4,y4] = rot_point(x0,y0,(180/p + th_m)*pi/180);

% archetto solo statore (pc + th_m)
mi_addarc(x1,y1,x2,y2,(pc + th_m),0.5);
[x,y] = rot_point(x0+eps,y0+eps,(-pc+th_m)/2*pi/180);
mi_selectarcsegment(x,y);
mi_setarcsegmentprop(0.5, 'APmove', 0, group);

% arco  comune (180/p - pc - th_m)
mi_addarc(x2,y2,x3,y3,(180/p-pc-th_m),0.5);
[x,y] = rot_point(x0+eps,y0+eps,(-pc+180/p+th_m)/2*pi/180);
mi_selectarcsegment(x,y);
mi_setarcsegmentprop(0.5, 'None', 0, group);

% archetto solo rotore (pc + th_m)
mi_addarc(x3,y3,x4,y4,(pc + th_m),0.5);
[x,y] = rot_point(x0+eps,y0+eps,(-pc+180/p+th_m+180/p)/2*pi/180);
mi_selectarcsegment(x,y);
mi_setarcsegmentprop(0.5, 'APmove', 0, group);


