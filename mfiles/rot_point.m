% rot_point - rotate the point (xi,yi) by the angle "angle"

function [xo,yo] = rot_point(xi,yi,angle)

xo = real((xi + 1j*yi).* exp(1j*angle));
yo = imag((xi + 1j*yi).* exp(1j*angle));
