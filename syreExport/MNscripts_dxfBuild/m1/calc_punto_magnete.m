% calc_punto_magnete.m
% calcola l'intersezione tra circonferenza centrata in x0,0 di raggio r1 e
% la retta che delimita il magnete, con fase gammam
% input: raggi r ed r1
% output: x,y del punto

function [x,y] = calc_punto_magnete(r1,gammam,x0)

x = x0 - r1 * cos(gammam);
y = r1 * sin(gammam);