%% modifiche 10 gen 2010 
% - salvo in geo le coordintate x_fe, y_fe di alcuni punti (punta dente,
% dente, giogo) e LDente (larg dente) per valutare le Pfe statore

master_seg_arcseg2
geo.wt = LDente;

% nodi da saltare (fanno sballare il resto e non servono a una mazza)
non_disegnare = [5 7 11 19  12 13 20];

for jj = 1:size(mast,1)
    if (sum(jj == non_disegnare) == 0)
        x1 = mast(jj,1); y1 = mast(jj,2);
        mi_addnode(x1,y1);
        mi_selectnode(x1,y1); mi_setnodeprop('None',1);
    end
%     keyboard
end

for jj = 1:size(arcseg,1)
    x1 = mast(arcseg(jj,1),1); y1 = mast(arcseg(jj,1),2);
    x2 = mast(arcseg(jj,2),1); y2 = mast(arcseg(jj,2),2);
    mi_drawarc(x1,y1,x2,y2,arcseg(jj,3),arcseg(jj,4));
    mi_selectarcsegment(mean([x1 x2]),mean([y1 y2]));
    mi_setarcsegmentprop(arcseg(jj,4), 'None', 0, group);
    % ultimo - aggiungo boundary
    if jj == size(arcseg,1)
        mi_setarcsegmentprop(arcseg(jj,4), 'A=0', 0, 1);
    end
%     keyboard
end

for jj = 1:size(segm,1)
    x1 = mast(segm(jj,1),1); y1 = mast(segm(jj,1),2);
    x2 = mast(segm(jj,2),1); y2 = mast(segm(jj,2),2);
    mi_drawline(x1,y1,x2,y2);
    mi_selectsegment(mean([x1 x2]),mean([y1 y2]));
    mi_setsegmentprop('None', res, 0, 0, group);
%     keyboard
end

mi_selectgroup(1);
mi_mirror2(0,0,10,0,4);
%% assegna materiali
% aria
mi_addblocklabel(mast(5,1),0);
mi_selectlabel(mast(5,1),0);
mi_setblockprop('Air', 0, res, 'None', 0, group, 1);
mi_clearselected;
% lamierino
% testa dente pos
x1 = mast(4,1); y1 = mast(4,2);
x2 = mast(9,1); y2 = mast(9,2);
temp_x = mean([x1 x2]); temp_y = mean([y1 y2]); 
mi_addblocklabel(temp_x,temp_y);
mi_selectlabel(temp_x,temp_y);
geo.x_fe = [geo.x_fe temp_x];
geo.y_fe = [geo.y_fe temp_y];
% testa dente neg
mi_addblocklabel(temp_x,-temp_y);
mi_selectlabel(temp_x,-temp_y);
geo.x_fe = [geo.x_fe -temp_x];
geo.y_fe = [geo.y_fe -temp_y];
% salvo un punto lungo il dente (uno per parte)
angolo = angle(temp_x + 1i * temp_y);
modulo = xr + g + lt/2;
[temp_x,temp_y] = pol2cart(angolo,modulo);
geo.x_fe = [geo.x_fe temp_x];
geo.y_fe = [geo.y_fe temp_y];
geo.x_fe = [geo.x_fe -temp_x];
geo.y_fe = [geo.y_fe -temp_y];
% giogo
x1 = mast(15,1); y1 = mast(15,2);
x2 = mast(14,1); y2 = mast(14,2);
temp_x = mean([x1 x2]); temp_y = mean([y1 y2]); 
mi_addblocklabel(temp_x,temp_y);
mi_selectlabel(temp_x,temp_y);
geo.x_fe = [geo.x_fe temp_x];
geo.y_fe = [geo.y_fe temp_y];
mi_setblockprop(steel, 0, res, 'None', 0, group, 1);
mi_clearselected
% rame
x1 = mast(17,1); y1 = mast(17,2);
x2 = mast(9,1); y2 = mast(9,2);
mi_addblocklabel(mean([x1 x2]),mean([y1 y2]));
mi_selectlabel(mean([x1 x2]),mean([y1 y2]));
x1 = mast(18,1); y1 = mast(18,2);
x2 = mast(14,1); y2 = mast(14,2);
mi_addblocklabel(mean([x1 x2]),mean([y1 y2]));
mi_selectlabel(mean([x1 x2]),mean([y1 y2]));
mi_setblockprop(copper, 0, res, 'None', 0, group, 1);
mi_clearselected

mi_zoomnatural

