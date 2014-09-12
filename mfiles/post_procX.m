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

function sol=post_procX(FemmProblem,ansfile,geo,th_m,thjj,id,iq)
myfpproc = fpproc();
myfpproc.setFemmProblem(FemmProblem);
myfpproc.opendocument(ansfile);
ns = geo.ns;
p = geo.p;
ps = geo.ps;
degs = 180/p*ps;
xr = geo.xr;
gap = geo.g;
pc  = 360/(ns*p)/2;

% load phase flux linkages
f1 = fluxLink(myfpproc,'fase1',p,ps);
f2 = fluxLink(myfpproc,'fase2',p,ps);
f3 = fluxLink(myfpproc,'fase3',p,ps);

% evaluate torque
% T1 - from the innermost integration line (rot + gap/6)
x = xr + gap*1/6;
ang0 = th_m; 
ang1 = degs + th_m;
T1=torque(myfpproc,x,ang0,ang1,geo);

% T2 - from the outermost integration line (stat - gap/6)
x = xr + gap*5/6;
ang0 = -pc; ang1 = degs-pc;
T2=torque(myfpproc,x,ang0,ang1,geo);

% T3 - from an intermediate line (rot + gap/2)
x = xr + gap*1/2;
ang0 = -pc; ang1 = degs-pc;
T3=torque(myfpproc,x,ang0,ang1,geo);

% dq flux linkaged
fdq = abc2dq(f1,f2,f3,thjj*pi/180);

% string of SOL
sol = [thjj id iq fdq(1) fdq(2) mean([T1,T2,T3])];

function f1=fluxLink(myfpproc,phaseName,p,ps)
temp_outX = myfpproc.getcircuitprops(phaseName);
temp_outXn= myfpproc.getcircuitprops([phaseName 'n']);
temp_out = temp_outX - temp_outXn;
f1 = temp_out(3) * 2 * p/ps;

function T1=torque(myfpproc,x,ang0,ang1,geo)
p = geo.p;
ps = geo.ps;
degs = 180/p*ps;
[x1,y1] = rot_point(x,0,ang0*pi/180);
[x2,y2] = rot_point(x,0,ang1*pi/180);
c=[x1 y1;x2 y2];
c=bendcountour(c,degs,0.5);
myfpproc.addcontour(c(:,1),c(:,2));
T1 = myfpproc.lineintegral(4);
T1 = T1(1) * 2 * p/ps;
myfpproc.clearcontour();
