% Matteo 2013/02/12
% retta per un punto tangente ad una circonferenza data.
% Retta tangente ad una circonferenza passante per un punto noto, di
% coefficiente angolare qualunque non noto
% r xp yp: rappresentano rispettivamente raggio e centro della
% circonferenza
% x1 y1 è il punto per cui passa la retta
% In output: eq retta mx+d, pto di tangenza x2,y2 
% x1=9/4; y1=0; xp=1; yp=0; r=1;

function [m d x2 y2]=retta_per_1_pto_tg_circonferenza(x1,y1,r,xp,yp)

m1=(-(y1-yp)*(xp-x1)+sqrt((y1-yp)^2*(xp-x1)^2-((xp-x1)^2-r^2)*((yp-y1)^2-r^2)))/((xp-x1)^2-r^2);
m2=(-(y1-yp)*(xp-x1)-sqrt((y1-yp)^2*(xp-x1)^2-((xp-x1)^2-r^2)*((yp-y1)^2-r^2)))/((xp-x1)^2-r^2);
m=[m1,m2];
d=y1-m.*x1;
a=-m; b=1; c=m.*x1-y1;
x2=(b^2*xp-a.*b.*yp-a.*c)./(a.^2+b^2);
y2=(a.^2.*yp-a.*b*xp-b.*c)./(a.^2+b^2);