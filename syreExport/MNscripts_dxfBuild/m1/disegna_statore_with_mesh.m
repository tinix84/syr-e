% disegna_statore.m - 24 08 09 - GMP

%% modifica 10 gennaio 2010 - gp
% - porto fuori geo per salvare i punti x_fe, y_fe di statore

function geo = disegna_statore_with_mesh(geo,m,fem)

group = 1;

ns = geo.ns;
b = geo.b;
kt = geo.kt;
acs = geo.acs;
avv = geo.avv;
Nbob = geo.Nbob;
lt = geo.lt;
p = geo.p;
r = geo.r;
xr = geo.xr;
g = geo.g;
Kbuccia=geo.kbuccia;
D0=geo.D0;
RaccordoFC=geo.RaccordoFC;
steel = m.steel_stat;   %different stator steel compared to rotor used; improve this using real different seel;
copper = m.copper;
% mesh resolution
res=fem.res;
res_traf=fem.res_traf;

pc = 360/(ns*p)/2;  % mezzo passo cava

mi_addcircprop('fase1', 0, 1);
mi_addcircprop('fase2', 0, 1);
mi_addcircprop('fase3', 0, 1);
mi_addcircprop('fase1n', 0, 1);
mi_addcircprop('fase2n', 0, 1);
mi_addcircprop('fase3n', 0, 1);

disegna_cava_with_mesh;
replica_cava;
% gap statore - strato 3
x1 = mast(4,1); y1 = -mast(4,2);
[x2,y2] = rot_point(xr+2/3*g,0,-pc*pi/180);
mi_drawline(x1,y1,x2,y2);
mi_selectsegment(mean([x1 x2]),mean([y1 y2]));
mi_setsegmentprop('APg3', res_traf, 0, 0, group);
mi_selectnode(x1,y1); mi_setnodeprop('None',group);
mi_selectnode(x2,y2); mi_setnodeprop('None',group);

[x1,y1] = rot_point(x1,y1,180/p*pi/180);
[x2,y2] = rot_point(x2,y2,180/p*pi/180);
mi_drawline(x1,y1,x2,y2);
mi_selectsegment(mean([x1 x2]),mean([y1 y2]));
mi_setsegmentprop('APg3', res, 0, 0, group);
mi_selectnode(x1,y1); mi_setnodeprop('None',group);
mi_selectnode(x2,y2); mi_setnodeprop('None',group);

[x1,y1] = rot_point(xr+2/3*g,0,-pc*pi/180);
mi_drawarc(x1,y1,x2,y2,180/p,pc/5);
[x1,y1] = rot_point(xr+4.5/6*g,0,90/p*pi/180);
mi_selectarcsegment(x1,y1);
mi_setarcsegmentprop(res_traf, 'None', 0, group);
mi_selectnode(x1,y1); mi_setnodeprop('None',group);
mi_selectnode(x2,y2); mi_setnodeprop('None',group);

% gap statore - strato 2
[x1,y1] = rot_point(xr+2/3*g,0,-pc*pi/180);
[x2,y2] = rot_point(xr+1/3*g,0,-pc*pi/180);
mi_drawline(x1,y1,x2,y2);
mi_selectsegment(mean([x1 x2]),mean([y1 y2]));
mi_setsegmentprop('APg2', res_traf, 0, 0, group);
mi_selectnode(x1,y1); mi_setnodeprop('None',group);
mi_selectnode(x2,y2); mi_setnodeprop('None',group);

[x1,y1] = rot_point(x1,y1,180/p*pi/180);
[x2,y2] = rot_point(x2,y2,180/p*pi/180);
mi_drawline(x1,y1,x2,y2);
mi_selectsegment(mean([x1 x2]),mean([y1 y2]));
mi_setsegmentprop('APg2', res, 0, 0, group);
mi_selectnode(x1,y1); mi_setnodeprop('None',group);
mi_selectnode(x2,y2); mi_setnodeprop('None',group);

% blocchi - strato 3
[x1,y1] = rot_point(xr+5/6*g,0,90/p*pi/180);
mi_addblocklabel(x1,y1);
mi_selectlabel(x1,y1);
mi_setblockprop('Air', 0, res_traf, 'None', 0, group, 1);

% strato 2
[x1,y1] = rot_point(xr+3/6*g,0,90/p*pi/180);
mi_addblocklabel(x1,y1);
mi_selectlabel(x1,y1);
mi_setblockprop('Air', 0, res_traf, 'None', 0, group, 1);

mi_clearselected


%% avvolgimento
% sel strato interno
x1 = mast(17,1); y1 = mast(17,2);
x2 = mast(9,1); y2 = mast(9,2);
xin0 = mean([x1 x2]); yin0 = mean([y1 y2]);
% sel strato esterno
x1 = mast(18,1); y1 = mast(18,2);
x2 = mast(14,1); y2 = mast(14,2);
xout0 = mean([x1 x2]); yout0 = mean([y1 y2]);
for jj = 1:ns/2
    % strato interno
    if avv(1,jj) > 0
        fase = ['fase' num2str(avv(1,jj))];
    else
        temp = num2str(avv(1,jj));
        fase = ['fase' temp(2) 'n'];
    end
    [x,y] = rot_point(xin0,yin0,(jj-1) * pc*2*pi/180);
    mi_selectlabel(x,y);
    mi_setblockprop(copper, 0, res, fase, 0, group, Nbob);
    mi_clearselected
    % strato esterno
    if avv(2,jj) > 0
        fase = ['fase' num2str(avv(2,jj))];
    else
        temp = num2str(avv(2,jj));
        fase = ['fase' temp(2) 'n'];
    end
    [x,y] = rot_point(xout0,yout0,(jj-1) * pc*2*pi/180);
    mi_selectlabel(x,y);
    mi_setblockprop(copper, 0, res, fase, 0, group, Nbob);
    mi_clearselected
end
