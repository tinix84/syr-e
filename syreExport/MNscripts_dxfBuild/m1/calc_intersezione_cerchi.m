% calc_intersezione_cerchi.m
% calcola l'intersezione tra circonferenza centrata in 0,0 di raggio r e
% circonferenza centrata in x0,0 di raggio r1
% input: raggi r ed r1
% output: x,y del punto

function [x,y] = calc_intersezione_cerchi(r,r1,x0)

x = (r.^2 - r1.^2 + x0^2)/(2*x0);
gamma = acos(x ./ r);
y = r .* sin(gamma);