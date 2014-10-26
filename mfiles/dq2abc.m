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

% trasformazione dq -> abc

function i123 = dq2abc(id,iq,theta)

% matrice di trasformazione (2 -> 3)
T23 = [1      0
    -0.5   sqrt(3)/2
    -0.5  -sqrt(3)/2];

% dq -> alpha beta
iab = (id + j * iq) .* exp(j*theta);
% alpha beta -> 123
i123 = T23 * [real(iab);imag(iab)];
