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
    case 'MO_GA'
        gamma = gamma_in;
        nsim = geo.nsim_MOOA;
        xdeg = geo.delta_sim_MOOA;
    case 'singt'
        nsim = geo.nsim_singt; xdeg = geo.delta_sim_singt;gamma = gamma_in;
    case 'singm'
        nsim = geo.nsim_singt; xdeg = geo.delta_sim_singt;gamma = gamma_in;       
end

% pathname = geo.pathname;
pathname = cd;

th0 = geo.th0;
p   = geo.p;
xr  = geo.xr;
gap = geo.g;
ns  = geo.ns;
pc  = 360/(ns*p)/2;
ps  = geo.ps;
% l = geo.l;
%% simulation angle
gradi_da_sim=180/p*ps;

id = io * cos(gamma * pi/180);
iq = io * sin(gamma * pi/180);

Hc = geo.Hc;

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

% for (single thread) or parfor (multi-thread) simulation cycle

% % if isOpen > 0
% %     
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     %%%%%%%%%%%%%%%%% PARFOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     
% %     % temporary files
% %     for i=1:32
% %         names_o{i}=['run' num2str(i)];
% %     end
% % keyboard   
% %     for jj = 1:nsim
% %         
% %         tmp_fem=[pathname,names_o{jj} '.fem'];
% % %         copyfile([pathname,'mot0.fem'],tmp_fem);
% %         h_temp=actxserver('femm.ActiveFEMM');
% %         callfemm_parfor([ 'setcurrentdirectory(' , quote(pathname) , ')'],h_temp);
% %         opendocument([pathname,'\mot0.femm']);
% %         
% %         th_m = (th(jj) - th0)/p;
% %         
% %         % assign the phase current values to the FEMM circuits
% %         i1 = i1_tmp(jj);
% %         i2 = i2_tmp(jj);
% %         i3 = i3_tmp(jj);
% %         mi_modifycircprop('fase1',1,i1);
% %         mi_modifycircprop('fase1n',1,-i1);
% %         mi_modifycircprop('fase2',1,i2);
% %         mi_modifycircprop('fase2n',1,-i2);
% %         mi_modifycircprop('fase3',1,i3);
% %         mi_modifycircprop('fase3n',1,-i3);
% %         % assign the Hc property to the bonded magnets
% %         mi_modifymaterial('Bonded-Magnet',3,Hc);
% %         % delete the airgap arc prior to moving the rotor
% %         mi_selectgroup(20), mi_deleteselectedarcsegments;
% %         % rotate the rotor
% %         mi_selectgroup(2), mi_moverotate(0,0,th_m);
% %         % redraw the airgap arc
% %         draw_airgap_arc_with_mesh(geo,th_m,geo.fem)
% %         
% %         mi_analyze(1);
% %         
% %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         %%%%%%%%%%%%%%%%% POST_PROC %%%%%%%%%%%%%%%%%%%%%%%%
% %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         
% %         mi_loadsolution_parfor(h_temp);
% %         
% %         % load phase flux linkages
% %         temp_out = mo_getcircuitproperties('fase1');
% %         temp_out = temp_out - mo_getcircuitproperties('fase1n');
% %         f1 = temp_out(3) * 2 * p/ps;
% %         temp_out = mo_getcircuitproperties('fase2');
% %         temp_out = temp_out - mo_getcircuitproperties('fase2n');
% %         f2 = temp_out(3) * 2 * p/ps;
% %         temp_out = mo_getcircuitproperties('fase3');
% %         temp_out = temp_out - mo_getcircuitproperties('fase3n');
% %         f3 = temp_out(3) * 2 * p/ps;
% %         
% %         % evaluate torque
% %         % T1 - from the innermost integration line (rot + gap/6)
% %         x = xr + gap*1/6;
% %         ang0 = th_m; ang1 = gradi_da_sim + th_m;
% %         [x1,y1] = rot_point(x,0,ang0*pi/180);
% %         [x2,y2] = rot_point(x,0,ang1*pi/180);
% %         mo_addcontour(x1,y1);
% %         mo_addcontour(x2,y2);
% %         mo_bendcontour(gradi_da_sim,0.5);
% %         
% %         T1 = mo_lineintegral(4);
% %         T1 = T1(1) * 2 * p/ps;
% %         mo_clearcontour();
% %         
% %         % T2 - from the outermost integration line (stat - gap/6)
% %         x = xr + gap*5/6;
% %         ang0 = -pc; ang1 = gradi_da_sim-pc;
% %         
% %         [x1,y1] = rot_point(x,0,ang0*pi/180);
% %         [x2,y2] = rot_point(x,0,ang1*pi/180);
% %         mo_addcontour(x1,y1);
% %         mo_addcontour(x2,y2);
% %         mo_bendcontour(gradi_da_sim,0.5);
% %         T2 = mo_lineintegral(4);
% %         T2 = T2(1) * 2 * p/ps;
% %         mo_clearcontour();
% %         
% %         % T3 - from an intermediate line (rot + gap/2)
% %         x = xr + gap*1/2;
% %         ang0 = -pc; ang1 = gradi_da_sim-pc;
% %         
% %         [x1,y1] = rot_point(x,0,ang0*pi/180);
% %         [x2,y2] = rot_point(x,0,ang1*pi/180);
% %         mo_addcontour(x1,y1);
% %         mo_addcontour(x2,y2);
% %         mo_bendcontour(gradi_da_sim,0.5);
% %         T3 = mo_lineintegral(4);
% %         T3 = T3(1) * 2 * p/ps;
% %         mo_clearcontour();
% %         
% %         % dq flux linkaged
% %         fdq = abc2dq(f1,f2,f3,th(jj)*pi/180);
% %         
% %         % string of SOL
% %         sol = [th(jj) id iq fdq(1) fdq(2) mean([T1,T2,T3])];
% %         
% %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         %%%%%%%%%%%%%%%%%% END OF POST PROC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         
% %         SOL_PARFOR_temp(jj,:) = sol;
% %         h_temp.delete;
% %         delete(tmp_fem);
% %         ans_temp=strrep(tmp_fem,'.fem','.ans');
% %         delete(ans_temp);
% %     end
% %     
% %     SOL = SOL_PARFOR_temp;
% %     
% % else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% ciclo for %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    for jj = 1:nsim
        
        th_m = (th(jj) - th0)/p;
       
        openfemm
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
        
        % assign the Hc property to the bonded magnets
        mi_modifymaterial('Bonded-Magnet',3,Hc);
        
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
