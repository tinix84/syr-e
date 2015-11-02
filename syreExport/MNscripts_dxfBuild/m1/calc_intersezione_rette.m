% calc_intersezione_rette
% calcola l'interersezione tra la retta1 passante per i punti x1,y1 x2,y2
% e la retta2 passante per i punti x3,y3 x4,y4

function [x,y] = calc_intersezione_rette(retta1,retta2)
x1 = retta1(1); y1 = retta1(2);
x2 = retta1(3); y2 = retta1(4);
x3 = retta2(1); y3 = retta2(2);
x4 = retta2(3); y4 = retta2(4);

% a-b: coefficienti della retta passante per 1-2
a = (y2 - y1) /(x2 - x1);
b = y1 - a * x1;
% c-d: coefficienti della retta passante per 3-4
c = (y4 - y3) /(x4 - x3);
d = y3 - c * x3;

% intersezione di y = ax +b e y = cx + d
x = (d - b) / (a - c);
y = a * x + b;

