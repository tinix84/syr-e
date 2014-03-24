%Calcolo area mediante prodotto vettoriale:
function [A]=calc_area(x1,x2,x3,y1,y2,y3)
vx=x2-x1; vy=y2-y1;
wx=x3-x1; wy=y3-y1;
A=abs(1/2*det([vx vy;wx wy]));