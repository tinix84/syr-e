%calcola intersezione tra 2 rette
%rette espresse in forma cartesiana.
function [x y]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2)
A=[a1 b1;a2 b2];
B=[-c1;-c2];
sol=A^-1*B;
x=sol(1,1);
y=sol(2,1);