% trasformazione dq -> abc

function i123 = dq2abc(id,iq,theta)

% matrice di trasformazione (2 -> 3)
T23 = [1      0
    -0.5   sqrt(3)/2
    -0.5  -sqrt(3)/2];

% dq -> alpha beta
iab = (id + j * iq) * exp(j*theta);
% alpha beta -> 123
i123 = T23 * [real(iab);imag(iab)];
