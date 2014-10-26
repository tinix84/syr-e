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

%% Matteo 31/10/2011
%Intersezione retta circonferenza:
function [x y]=intersezione_retta_circonferenza(x0,y0,r)

[v]=roots([1, 1*(y0-x0), 1/2*(y0^2+x0^2-2*y0*x0-r^2)]);

t=1;
for i=1:length(v)
    if v(i)>=0
        w(t)=v(i);
    else
        
    end
    t=t+1;
end
x=w;
y=y0+x-x0;

