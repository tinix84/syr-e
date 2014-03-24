% replica cava

mi_selectgroup(1);
mi_copyrotate2(0, 0, s.AS1*2, ns/2-1, 4 );
mi_clearselected;


% anti periodicità
% espansione dente
x1 = mast(4,1); y1 = -mast(4,2);
x2 = mast(10,1); y2 = -mast(10,2);
mi_selectsegment(mean([x1 x2]),mean([y1 y2]));
[x1, y1] = rot_point(x1,y1,s.AS1*2*ns/2*pi/180);
[x2, y2] = rot_point(x2,y2,s.AS1*2*ns/2*pi/180);
mi_selectsegment(mean([x1 x2]),mean([y1 y2]));
mi_setsegmentprop('APs1', 0, 1, 0, group);
mi_clearselected;

% espansione dente
x1 = mast(16,1); y1 = -mast(16,2);
x2 = mast(10,1); y2 = -mast(10,2);
mi_selectsegment(mean([x1 x2]),mean([y1 y2]));
[x1, y1] = rot_point(x1,y1,s.AS1*2*ns/2*pi/180);
[x2, y2] = rot_point(x2,y2,s.AS1*2*ns/2*pi/180);
mi_selectsegment(mean([x1 x2]),mean([y1 y2]));
mi_setsegmentprop('APs2', 0, 1, 0, group);
mi_clearselected;

