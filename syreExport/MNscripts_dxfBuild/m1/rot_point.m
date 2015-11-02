% rot_point - rotate the point (xi,yi) by the angle "angle"

function [xo,yo] = rot(xi,yi,angle)

xo = real((xi + j*yi).* exp(j*angle));
yo = imag((xi + j*yi).* exp(j*angle));
