% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

% abc -> dq general transformation for multi 3-phase machines
% same results as dq2abc.m for the case n3phase=1

function idq = abc2dq(i1,i2,i3,theta,n3phase,index)

% index vary between 0 to (n3phase-1)
% n3phase = number of 3-phase circuits

if (nargin == 4)
    theta_i = 0;
else
    delta = pi/3;                   %sextant amplitude
    theta_i = delta*index/n3phase;  %angular difference between alpha-axis and first phase of the n3phase-th 3-phase circuit
end    
% matrice di trasformazione (3 -> 2)
T32 = 2/3 * [cos(theta_i)     cos(theta_i+2*pi/3)    cos(theta_i-2*pi/3)
    sin(theta_i)     sin(theta_i+2*pi/3)    sin(theta_i-2*pi/3)];

% 123 -> alpha beta
iab = T32 * [i1;i2;i3];
% dq -> alpha beta
temp = (iab(1,:) + j * iab(2,:)) .* exp(-j*theta);

idq = [real(temp);imag(temp)];


