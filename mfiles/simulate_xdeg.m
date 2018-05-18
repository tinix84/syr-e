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
% 
% Update 02/2018: SOL is a structure.
%function [SOL] = simulate_xdeg(geo,io,Br,gamma_in,eval_type)
function [SOL] = simulate_xdeg(geo,per,eval_type,pathname,filename)
% NB: pathname with final slash
% number of simulation that must be done respect to eval type
switch eval_type
    case 'MO_OA'
        gamma = per.gamma;
        nsim = geo.nsim_MOOA;
        xdeg = geo.delta_sim_MOOA;
        randFactor = geo.randFactor;
    case 'MO_GA'
        gamma = per.gamma;
        nsim = geo.nsim_MOOA;
        xdeg = geo.delta_sim_MOOA;
        randFactor = geo.randFactor;
    case 'singt'
        nsim = geo.nsim_singt;
        xdeg = geo.delta_sim_singt;
        gamma = per.gamma;
        randFactor = 0;
    case 'singm'
        nsim = geo.nsim_singt;
        xdeg = geo.delta_sim_singt;
        gamma = per.gamma;
        randFactor = 0;
end

%pathname = pwd();

th0 = geo.th0;
p   = geo.p;
r  = geo.r;
gap = geo.g;
ns  = geo.ns;
pc  = 360/(ns*p)/2;
ps  = geo.ps;
n3phase = geo.n3phase; %AS number of 3-phase circuits
% simulation angle
gradi_da_sim=180/p*ps;

Hc = per.BrPP/(4e-7*pi);

%SOL = [];

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
iAmp = per.overload*calc_io(geo,per);
iAmpCoil = iAmp*geo.Nbob*geo.n3phase; %AS

id = iAmpCoil * cos(gamma * pi/180);
iq = iAmpCoil * sin(gamma * pi/180);
i_tmp = zeros(3*n3phase,nsim);   %matrix containing all phase current values for the simulated rotor position
for ik=0:(n3phase-1)
    for ij=1:nsim
        if geo.avv_flag((3*ik)+1)==1 && geo.avv_flag((3*ik)+2)==1 && geo.avv_flag((3*ik)+3)==1
            i123 = dq2abc(id,iq,th(ij)*pi/180,n3phase,ik);
            i_tmp((3*ik)+1,ij) = i123(1);
            i_tmp((3*ik)+2,ij) = i123(2);
            i_tmp((3*ik)+3,ij) = i123(3);
        else
            if geo.avv_flag((3*ik)+1)==0 && geo.avv_flag((3*ik)+2)==0 && geo.avv_flag((3*ik)+3)==0
                i_tmp((3*ik)+1,ij) = 0;
                i_tmp((3*ik)+2,ij) = 0;
                i_tmp((3*ik)+3,ij) = 0;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% ciclo for %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Open and draw motor once, rotate and simulate nsim positions by Chao Lu 17/01/2017
openfemm
%main_minimize
opendocument([pathname,filename]);

SOL.th = zeros(1,nsim); % electrical angle in degree
SOL.id = zeros(1,nsim); % d-axis current
SOL.iq = zeros(1,nsim); % q-axis current
SOL.fd = zeros(1,nsim); % d-axis flux linkage
SOL.fq = zeros(1,nsim); % q-axis flux linkage
SOL.T  = zeros(1,nsim); % get from the line integral along the airgap (standard)
%SOL.Tb = zeros(1,nsim); % get from block integral (suggested by David Meeker)
SOL.F  = zeros(1,nsim); % radial and tangential force of the simulated portion

for jj = 1:nsim

    th_m = (th(jj) - th0)/p;
    
    %         openfemm
    %         %main_minimize
    %         opendocument([pathname,'\mot0.fem']);
    
    % assign the phase current values to the FEMM circuits
    for ik=0:(n3phase-1)
        phase_name{3*ik+1}=strcat('fase',num2str(3*ik+1));
        phase_name{3*ik+2}=strcat('fase',num2str(3*ik+2));
        phase_name{3*ik+3}=strcat('fase',num2str(3*ik+3));
        phase_name_neg{3*ik+1}=strcat('fase',num2str(3*ik+1),'n');
        phase_name_neg{3*ik+2}=strcat('fase',num2str(3*ik+2),'n');
        phase_name_neg{3*ik+3}=strcat('fase',num2str(3*ik+3),'n');
        
        % change current value in FEMM
        mi_modifycircprop(phase_name{3*ik+1}, 1,i_tmp((3*ik)+1,jj));
        mi_modifycircprop(phase_name{3*ik+2}, 1,i_tmp((3*ik)+2,jj));
        mi_modifycircprop(phase_name{3*ik+3}, 1,i_tmp((3*ik)+3,jj));
        mi_modifycircprop(phase_name_neg{3*ik+1}, 1,-i_tmp((3*ik)+1,jj));
        mi_modifycircprop(phase_name_neg{3*ik+2}, 1,-i_tmp((3*ik)+2,jj));
        mi_modifycircprop(phase_name_neg{3*ik+3}, 1,-i_tmp((3*ik)+3,jj));
    end
    
    % assign the Hc property to each of the bonded magnets
    if strcmp(geo.RotType,'SPM')
        mi_modifymaterial('Bonded-Magnet',3,Hc);
    else
        if length(Hc)==1
            Hc_vect = Hc*ones(1,length(geo.BLKLABELS.rotore.BarName));
        else
            Hc_vect=[Hc Hc];
        end
        for ii = 1:length(Hc_vect)
            mi_modifymaterial(['Bonded-Magnet_' num2str(ii)],3,Hc_vect(ii));
        end
    end
    
    % delete the airgap arc prior to moving the rotor
    mi_selectgroup(20), mi_deleteselectedarcsegments;
    % rotate the rotor
    mi_selectgroup(22),  mi_selectgroup(2)
    if jj > 1
        mi_moverotate(0,0,(th(jj) - th(jj-1))/p);
    else
        mi_moverotate(0,0,th_m);
    end
    % redraw the airgap arc
    if (ps<2*p)
        draw_airgap_arc_with_mesh(geo,th_m,geo.fem)
    else
        draw_airgap_arc_with_mesh_fullMachine(geo,th_m,geo.fem)
    end
    mi_saveas([pathname,'\mot_temp.fem']);
    mi_analyze(1);
    mi_loadsolution;
    post_proc;

    
    SOL.th(jj) = th(jj);
    SOL.id(jj) = id;
    SOL.iq(jj) = iq;
    SOL.fd(jj) = fd;  %AS
    SOL.fq(jj) = fq;  %AS
    SOL.T(jj)  = Tblock;
    %SOL.T(jj)  = mean([T1,T2,T3]);
    %SOL.Tb(jj) = Tblock;
    %SOL.F(jj)  = Frt;
% Added in struct VolPM - rev.Gallo
    SOL.VolPM(jj)= VolPM;

    %SOL = [SOL; sol];
    
end

mo_close, mi_close
closefemm

% Effective value of flux and current, simulation are done with one turns
% in slot and consequently, current in fem simulation is increase by the number of conductors in slot Nbob....

SOL.id=SOL.id/geo.Nbob/n3phase;
SOL.iq=SOL.iq/geo.Nbob/n3phase;
SOL.fd=SOL.fd*geo.Nbob*n3phase;
SOL.fq=SOL.fq*geo.Nbob*n3phase;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       FINE VECCHIO CICLO FOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % end
%% VERIFICA CHE I DUE CICLI FOR E PARFOR DIANO LO STESSO RISULTATO %%
% ERROR = SOL-SOL_PARFOR_temp
