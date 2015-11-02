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


% Input: nsim, id, iq;
% Output:
% - vector of the solutions SOL (nsim x 6);
% - Temp\sim_gamma_numerosimulazione.fem;
% - Temp\sim_gamma.mat (memorizza SOL)

% The nsim x 6 matrix is organized by columns as folows:
% SOL(1,:) = gamma,
% SOL(2,:) = id,
% SOL(3,:) = iq,
% SOL(4,:) = fd,
% SOL(5,:) = fq,
% SOL(6,:) = torque;

function [SOL] = simulate_xdeg(geo,io,gamma_in,eval_type)

% number of simulation that must be done respect to eval type
switch eval_type
    case 'MO_OA'
        gamma = gamma_in;
        nsim = geo.nsim_MOOA;
        xdeg = geo.delta_sim_MOOA;
        randFactor = geo.randFactor;
    case 'MO_GA'
        gamma = gamma_in;
        nsim = geo.nsim_MOOA;
        xdeg = geo.delta_sim_MOOA;
        randFactor = geo.randFactor;
    case 'singt'
        nsim = geo.nsim_singt; xdeg = geo.delta_sim_singt;gamma = gamma_in;
        randFactor = 0;
    case 'singm'
        nsim = geo.nsim_singt; xdeg = geo.delta_sim_singt;gamma = gamma_in;
        randFactor = 0;
end

pathname = cd;

th0 = geo.th0;
p   = geo.p;
r  = geo.r;
gap = geo.g;
ns  = geo.ns;
pc  = 360/(ns*p)/2;
ps  = geo.ps;

% simulation angle
gradi_da_sim=180/p*ps;

id = io * cos(gamma * pi/180);
iq = io * sin(gamma * pi/180);

% Hc = geo.Hc;
Hc = geo.Br/(4e-7*pi);

SOL = [];

% rotor positions
if strcmp(eval_type,'MO_OA')||strcmp(eval_type,'MO_GA')
    % during optimization, random position offset
    sim_step=xdeg/(nsim+0.5);
    offset=1*sim_step*rand;
    isOpen=0; %Disable parFor during optimization
else
    % during re-evaluation, regular position steps
    sim_step=xdeg/(nsim);
    offset=0;
    isOpen=1; %Enable parFor during optimization
end

teta=offset:sim_step:xdeg+offset;
teta=teta+(0.5-rand(1,length(teta)))*randFactor;


% disregard the last position
th=th0+[teta(1:nsim) teta(1)];

% evaluation of the phase current values for all positions to be simulated
i1_tmp = zeros(1,nsim); i2_tmp = i1_tmp; i3_tmp = i1_tmp;
for ij=1:nsim
    i123 = dq2abc(id,iq,th(ij)*pi/180);
    i1_tmp(ij) = i123(1);
    i2_tmp(ij) = i123(2);
    i3_tmp(ij) = i123(3);
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% ciclo for %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    for jj = 1:nsim
        
        th_m = (th(jj) - th0)/p;
       
        openfemm
        %main_minimize
        opendocument([pathname,'\mot0.fem']);
        
        % assign the phase current values to the FEMM circuits
        i1 = i1_tmp(jj);
        i2 = i2_tmp(jj);
        i3 = i3_tmp(jj);
        mi_modifycircprop('fase1',1,i1);
        mi_modifycircprop('fase1n',1,-i1);
        mi_modifycircprop('fase2',1,i2);
        mi_modifycircprop('fase2n',1,-i2);
        mi_modifycircprop('fase3',1,i3);
        mi_modifycircprop('fase3n',1,-i3);
        
        % assign the Hc property to each of the bonded magnets
        for j = 1:length(Hc)
            mi_modifymaterial(['Bonded-Magnet' num2str(j)],3,Hc(j));
        end
        
        % delete the airgap arc prior to moving the rotor
        mi_selectgroup(20), mi_deleteselectedarcsegments;
        
        % rotate the rotor
        mi_selectgroup(2), mi_moverotate(0,0,th_m);
        
        % redraw the airgap arc
        draw_airgap_arc_with_mesh(geo,th_m,geo.fem)
        
        mi_saveas([pathname,'\mot_temp.fem']);    
        mi_analyze(1);
        mi_loadsolution;
        post_proc;
        mo_close, mi_close
        closefemm
        
        SOL = [SOL; sol];
		
    end
    
    % Effective value of flux and current, simulation are done with one turns
    % in slot and consequently, current in fem simulation is increase by the number of conductors in slot Nbob....
    SOL(:,2)=SOL(:,2)/geo.Nbob;
    SOL(:,3)=SOL(:,3)/geo.Nbob;
    SOL(:,4)=SOL(:,4)*geo.Nbob;
    SOL(:,5)=SOL(:,5)*geo.Nbob;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       FINE VECCHIO CICLO FOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % end
%% VERIFICA CHE I DUE CICLI FOR E PARFOR DIANO LO STESSO RISULTATO %%
% ERROR = SOL-SOL_PARFOR_temp
