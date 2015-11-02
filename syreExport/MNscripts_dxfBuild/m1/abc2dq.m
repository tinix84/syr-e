% trasformazione abc -> dq
% i1,i2,i3 vettori riga


function idq = abc2dq(i1,i2,i3,theta)

% matrice di trasformazione (3 -> 2)    
T32 = 2/3 * [	1 	-0.5 		-0.5
                0    sqrt(3)/2    -sqrt(3)/2];

% 123 -> alpha beta
iab = T32 * [i1;i2;i3];
% dq -> alpha beta
temp = (iab(1) + j * iab(2)) * exp(-j*theta);

idq = [real(temp) imag(temp)];
