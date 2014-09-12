% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function [MatrWin] = MatrixWin(Q,p,isSingle,yq)
%% ===== MatrixWin ========================================================
%%              INPUT: Q = number of slots
%%                     p = number of pole pairs   
%%                     isSingle = flag for single layer winding
%%                     yq = coil throw
%%             OUTPUT: MatrWin = Matrix of windings (2xQ/periodicity) 
%% ========================================================================

m = 3; % number of phases
t = gcd(Q,p); % machine periodicity and max. number of parallel circuits
alfase = (2*pi/Q)*p; % electrical angle between slot
c = Q/m*t - fix(Q/m*t); % check if the ratio is entire
if c ~= 0
    disp('Winding not feasible');
    return
end
fase = zeros(Q,1); % matrix of the phases
d = 0;
for h = 1 : 1 : Q
    while d > 2*pi
        d = d - 2*pi;
    end
    if abs(d) < 1e-2 || abs(d - 2*pi) < 1e-2
        d = 0;
    end
    fase(h,:) = d; 
    d = d + alfase; % angles
end
slot = 1 : 1 : Q; % matrix of the slots
Y = diag(slot,0)*exp(1i*fase); % matrix of the spokes and of the phases
settore = pi/m; % angular sector 
A = zeros(m,Q); % matrix of the first layer (double-layer winding)
delta_set = 0.3*pi/180;
for h = 1 : 1 : m
    y = -settore/2 + (h-1)*2*settore + delta_set;
    x = settore/2 + (h-1)*2*settore + delta_set;
    z = y + pi;
    k = x + pi;
    if y >= 2*pi
        y = y - 2*pi;
        x = x - 2*pi;
    elseif z >= 2*pi
        z = z - 2*pi;
        k = k - 2*pi;
    end
    for j = 1 : 1 : Q
        if angle(Y(j)) >= 0
            d = angle(Y(j));
        elseif angle(Y(j)) < 0
            d = angle(Y(j)) + 2*pi;
        end
        if h == 1 && d >= (y + 2*pi - delta_set)
            d = d - 2*pi;
        end        
        if d > (y) && d < (x)
            A(h,j) = 1*h;
        elseif d > (z) && d < (k)
            A(h,j) = -1*h;
        end
    end
end
% A build choosing as reference the left side of the coils;
B = zeros(m,Q);
for h = 1 : 1 : m
    for j = 1 : 1 : Q
        if yq > 0
            c = j + yq;
            if c > Q
                c = (j + yq) - Q;
            end
            B(h,c) = - A(h,j);
        elseif yq < 0
            c = j + yq;
            if c < 1
                c = Q - abs(j + yq);
            elseif c == 0
                c = Q;
            end
            B(h,c) = -A(h,j);            
        end
    end
end
% From double-Layer to single-Layer
E = A;
% in the single-layer winding, one delets the odd spokes 
for i = 1 : m
    for j = 1 : Q
        if mod(j,2) == 0
            E(i,j) = 0;
        end
    end
end
F = zeros(m,Q);
for h = 1 : 1 : m
    for j = 1 : 1 : Q
        if yq > 0
            c = j + yq;
            if c > Q
                c = (j + yq) - Q;
            end
            F(h,c) = -E(h,j);
        elseif yq < 0
            c = j + yq;
            if c < 1
                c = Q - abs(j + yq);
            elseif c == 0
                c = Q;
            end
            F(h,c) = -E(h,j);            
        end
    end
end

I = E + F; % winding matrix in single-layer case

if isSingle == 0 %(double-layer case)
    
    MatrWin(1,:) = sum(A,1); % matrix of the first layer 
    
    MatrWin(2,:) = sum(B,1); % matrix of the second layer
    
    
else %(single-layer case)
    
    MatrWin(1,:) = sum(I,1); % matrix of the first layer
    
    MatrWin(2,:) = MatrWin(1,:); % matrix of the second layer
    
end
    
MatrWin = MatrWin(:,1:Q/t);
end

