% calc_distanza_punti.m
% calcola la distanza tra due punti di coordinate note

function d = calc_distanza_punti(punto1,punto2)

x1 = punto1(1); y1 = punto1(2);
x2 = punto2(1); y2 = punto2(2);

d = sqrt((x1 - x2)^2 + (y1 - y2)^2);