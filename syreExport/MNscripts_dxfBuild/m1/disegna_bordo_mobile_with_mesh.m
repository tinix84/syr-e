%% Disegna la linea di confine statore - rotore che va rifatta ad ogni posizione

function disegna_bordo_mobile_with_mesh(geo,th_m,fem)

global tipo_valutazione

mi_clearselected

p  = geo.p;
xr = geo.xr;
g  = geo.g;

ns = geo.ns;
% Mezzo passo cava
pc = 360/(ns*p)/2;

res=fem.res_traf;

gradi_da_sim = 180/p;  % si fa riferimento ai gradi meccanici della porzione di rot simulata

group = 20;

x0 = xr+1/3*g;
y0 = 0;

if (strcmp(tipo_valutazione,'gambj'))
    res_deg=res;     %NON è un vero res_deg anc se mantiene il nome
else
    res_deg=res*180/pi/x0;
%     res_deg=0.5;
end

angoli_bordo_mobile = [-pc, th_m, (-pc + gradi_da_sim), (gradi_da_sim+th_m)];
angoli_bordo_mobile = sort(angoli_bordo_mobile);
delta_angoli_bordo_mobile = diff(angoli_bordo_mobile);

[x1,y1] = rot_point(x0,y0,angoli_bordo_mobile(1)*pi/180);
[x2,y2] = rot_point(x0,y0,angoli_bordo_mobile(2)*pi/180);
[x3,y3] = rot_point(x0,y0,angoli_bordo_mobile(3)*pi/180);
[x4,y4] = rot_point(x0,y0,angoli_bordo_mobile(4)*pi/180);

num_seg=delta_angoli_bordo_mobile(1)*pi/180*x0/res;
mi_addarc(x1,y1,x2,y2,delta_angoli_bordo_mobile(1),num_seg);
[x,y] = rot_point(x1,y1,0.5*delta_angoli_bordo_mobile(1)*pi/180);
mi_selectarcsegment(x,y);
mi_setarcsegmentprop(res_deg, 'APmove', 0, group);

num_seg=delta_angoli_bordo_mobile(2)*pi/180*x0/res;
mi_addarc(x2,y2,x3,y3,delta_angoli_bordo_mobile(2),num_seg);
[x,y] = rot_point(x2,y2,0.5*delta_angoli_bordo_mobile(2)*pi/180);
mi_selectarcsegment(x,y);
mi_setarcsegmentprop(res_deg, 'None', 0, group);
mi_clearselected

num_seg=delta_angoli_bordo_mobile(3)*pi/180*x0/res;
mi_addarc(x3,y3,x4,y4,delta_angoli_bordo_mobile(3),num_seg);
[x,y] = rot_point(x3,y3,0.5*delta_angoli_bordo_mobile(3)*pi/180);
mi_selectarcsegment(x,y);
mi_setarcsegmentprop(res_deg, 'APmove', 0, group);
mi_clearselected
end



