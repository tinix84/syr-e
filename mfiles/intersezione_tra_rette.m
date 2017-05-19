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

%calcola intersezione tra 2 rette
%rette espresse in forma cartesiana.
function [x,y]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2)
A=[a1 b1;a2 b2];
B=[-c1;-c2];
sol=A^-1*B;
x=sol(1,1);
y=sol(2,1);


%% ax + by + c = 0
% x = -c1./a1 - b1*(a1*c2-c1*a2)./(a1*b1*a2-b2*a1^2); 
% y = (a1.*c2-c1.*a2)./(b1.*a2-a1.*b2);

