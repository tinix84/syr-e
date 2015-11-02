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


l = geo.l;
p = geo.p;

% load phase flux linkages
temp_out = mo_getcircuitproperties('fase1');
temp_out = temp_out - mo_getcircuitproperties('fase1n');
f1 = temp_out(3) * 2 * p/ps;
temp_out = mo_getcircuitproperties('fase2');
temp_out = temp_out - mo_getcircuitproperties('fase2n');
f2 = temp_out(3) * 2 * p/ps;
temp_out = mo_getcircuitproperties('fase3');
temp_out = temp_out - mo_getcircuitproperties('fase3n');
f3 = temp_out(3) * 2 * p/ps;

% evaluate torque
% T1 - from the innermost integration line (rot + gap/6)
x = r + gap*1/6;
ang0 = th_m; ang1 = gradi_da_sim + th_m;
[x1,y1] = rot_point(x,0,ang0*pi/180);
[x2,y2] = rot_point(x,0,ang1*pi/180);
mo_addcontour(x1,y1);
mo_addcontour(x2,y2);
mo_bendcontour(gradi_da_sim,0.5);

T1 = mo_lineintegral(4);
T1 = T1(1) * 2 * p/ps;
mo_clearcontour();

% T2 - from the outermost integration line (stat - gap/6)
x = r + gap*5/6;
ang0 = -pc; ang1 = gradi_da_sim-pc;

[x1,y1] = rot_point(x,0,ang0*pi/180);
[x2,y2] = rot_point(x,0,ang1*pi/180);
mo_addcontour(x1,y1);
mo_addcontour(x2,y2);
mo_bendcontour(gradi_da_sim,0.5);
T2 = mo_lineintegral(4);
T2 = T2(1) * 2 * p/ps;
mo_clearcontour();

% T3 - from an intermediate line (rot + gap/2)
x = r + gap*1/2;
ang0 = -pc; ang1 = gradi_da_sim-pc;

[x1,y1] = rot_point(x,0,ang0*pi/180);
[x2,y2] = rot_point(x,0,ang1*pi/180);
mo_addcontour(x1,y1);
mo_addcontour(x2,y2);
mo_bendcontour(gradi_da_sim,0.5);
T3 = mo_lineintegral(4);
T3 = T3(1) * 2 * p/ps;
mo_clearcontour();

% dq flux linkaged
fdq = abc2dq(f1,f2,f3,th(jj)*pi/180);

% MAURIZIO --> Acquisizione dei valori di Bx e By puntuali per effettuare
% il calcolo delle perdite nel ferro

nSim = geo.nsim_singt;              % Number of simulations
deltasim = geo.delta_sim_singt; 
tetaMAU = nSim/deltasim;
thMAU = th_m;
indMAU = jj;

% Storage of flux densities in all the mesh triangle
% keyboard
if indMAU == 1
    nnM = mo_numelements;               % Number of mesh elements
    bM = zeros(nSim,nnM);               % Matrix that will hold the Bx By info
    zM = zeros(nnM,1);                 % Matrix that will hold the mesh elements centroid coordinates as complex number
    aM = zeros(nnM,1);                 % Matrix that will hold the mesh elements area
    gM = zeros(nnM,1);                 % Matrix that will hold the mesh elements group number
    for mM = 1:nnM
        elm = mo_getelement(mM);
        zM(mM) = elm(4)+j*elm(5);
        aM(mM) = elm(6);
        gM(mM) = elm(7);
    end
end

uMAU = exp(j*thMAU*pi/180);            % It is a parameter that considers the rotor angle

% % keyboard

for mM = 1:nnM
    if (gM(mM) == 2)
        prot = zM(mM)*uMAU;
        bM(indMAU,mM) = (mo_getb(real(prot),imag(prot))*[1;j])/uMAU;
    end
    if (gM(mM) == 1)
        prot = zM(mM);
        bM(indMAU,mM) = (mo_getb(real(prot),imag(prot))*[1;j]);
    end
end
% keyboard

% string of SOL

sol = [th(jj) id iq fdq(1) fdq(2) mean([T1,T2,T3])];
fluxdens = [zM aM gM bM'];      % Saving all the info necessary to calculate Iron Losses





