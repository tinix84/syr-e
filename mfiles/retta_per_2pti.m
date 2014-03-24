%restituisce i parametri a b c di una retta passante per 2 pti:
function [a b c]=retta_per_2pti(x1,y1,x2,y2)
a=y1-y2;
b=x2-x1;
c=x1*(y2-y1)-y1*(x2-x1);