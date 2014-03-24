% calc_apertura_cerchio.m
% calcola l'apertura angolare rispetto al centro (x0,0) di un punto con
% coordinate polari (r, alpha) rispetto al centro (0,0)
% input: raggi r, alpha, x0
% output: angolo rispetto a (x0,0)

function beta = calc_apertura_cerchio(alpha,r,x0)

beta = atan(r * sin(alpha) ./ (x0 - r * cos(alpha)));