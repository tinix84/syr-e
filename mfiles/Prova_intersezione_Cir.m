%% Intersezione tra due circonferenze di raggio e centro qualunque; le
%% circonferenze vengono espresse in funzione delle coordinate polari
%% centro e raggio e successivamente trasformate in forma implicita per la
%% risoluzione del problema.

% a1=0, b1=-2, c1=1-4,  a2=-4, b2=0, c2=+4-4,

function [xInt1,yInt1,xInt2,yInt2]=Prova_intersezione_Cir(x1,y1,r1,x2,y2,r2)

a1=-2*x1; b1=-2*y1; c1=x1^2+y1^2-r1^2;  a2=-2*x2; b2=-2*y2; c2=x2^2+y2^2-r2^2;

A=((b2-b1)^2+(a2-a1)^2);
B=(a2-a1)*(c2-c1)+a1/2*(b2-b1)^2-b1/2*(b2-b1)*(a2-a1);
C=(c2-c1)^2-b1*(b2-b1)*(c2-c1)+c1*(b2-b1)^2;
xInt1=1/A*(-B+sqrt(B^2-A*C));
xInt2=1/A*(-B-sqrt(B^2-A*C));
yInt1=(-(a2-a1)*xInt1-(c2-c1))/(b2-b1);
yInt2=(-(a2-a1)*xInt2-(c2-c1))/(b2-b1);