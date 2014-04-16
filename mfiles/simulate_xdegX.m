
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

function [SOL,FemmProblem] = simulate_xdegX(FemmProblem,geo,nsim,xdeg,io,gamma,eval_type,boundnameAPmove)

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
if strcmp(eval_type,'MO_OA')
    % during optimization, random position offset
    sim_step=xdeg/(nsim-1+0.5);
    offset=1*sim_step*rand;
    isOpen=0; %Disable parFor during optimization
else
    % during re-evaluation, regular position steps
    sim_step=xdeg/(nsim-1);
    offset=0;
    isOpen=1; %Enable parFor during optimization
end


teta=0:sim_step:xdeg+offset;

% disregard the last position
th=th0+[teta(1:nsim-1) teta(1)];

% evaluation of the phase current values for all positions to be simulated
i1_tmp = zeros(1,nsim-1); i2_tmp = i1_tmp; i3_tmp = i1_tmp;
for ij=1:nsim-1
    i123 = dq2abc(id,iq,th(ij)*pi/180);
    i1_tmp(ij) = i123(1);
    i2_tmp(ij) = i123(2);
    i3_tmp(ij) = i123(3);
end

% for (single thread) or parfor (multi-thread) simulation cycle

if isOpen > 0
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% PARFOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % temporary files
    for i=1:32
        names_o{i}=['run' num2str(i)];
    end
    
    parfor jj = 1:nsim-1
        
        tmp_fem=[pathname,'\',names_o{jj} '.fem'];
        copyfile([pathname,'\mot0.fem'],tmp_fem);
        h_temp=actxserver('femm.ActiveFEMM');
        callfemm_parfor([ 'setcurrentdirectory(' , quote(pathname) , ')'],h_temp);
        opendocument_parfor(tmp_fem,h_temp);
        
        th_m = (th(jj) - th0)/p;
        
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
        
        mi_analyze(1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%% POST_PROC %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        mi_loadsolution_parfor(h_temp);
        
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
        x = xr + gap*1/6;
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
        x = xr + gap*5/6;
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
        x = xr + gap*1/2;
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
        
        % string of SOL
        sol = [th(jj) id iq fdq(1) fdq(2) mean([T1,T2,T3])];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% END OF POST PROC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        SOL_PARFOR_temp(jj,:) = sol;
        h_temp.delete;
        delete(tmp_fem);
        ans_temp=strrep(tmp_fem,'.fem','.ans');
        delete(ans_temp);
    end
    
    SOL = SOL_PARFOR_temp;
    
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% ciclo for %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for jj = 1:nsim-1
        th_m = (th(jj) - th0)/p;
        
        %opendocument([pathname,'\mot0.fem']);
        
        % assign the phase current values to the FEMM circuits
        i1 = i1_tmp(jj);
        i2 = i2_tmp(jj);
        i3 = i3_tmp(jj);
        FemmProblem.Circuits(1).TotalAmps_re=i1;
        FemmProblem.Circuits(2).TotalAmps_re=i2;
        FemmProblem.Circuits(3).TotalAmps_re=i3;
        FemmProblem.Circuits(4).TotalAmps_re=-i1;
        FemmProblem.Circuits(5).TotalAmps_re=-i2;
        FemmProblem.Circuits(6).TotalAmps_re=-i3;

%         mi_modifycircprop('fase1',1,i1);
%         mi_modifycircprop('fase1n',1,-i1);
%         mi_modifycircprop('fase2',1,i2);
%         mi_modifycircprop('fase2n',1,-i2);
%         mi_modifycircprop('fase3',1,i3);
%         mi_modifycircprop('fase3n',1,-i3);
        
        % assign the Hc property to the bonded magnets
        %mi_modifymaterial('Bonded-Magnet',3,Hc);
        M=newmaterial_mfemm('Bonded-Magnet','H_c',Hc);
        FemmProblem = addmaterials_mfemm(FemmProblem,M);
        % delete the airgap arc prior to moving the rotor
        %mi_selectgroup(20), mi_deleteselectedarcsegments;
        FemmProblem = deletegroup_mfemm(FemmProblem, 20);
        % rotate the rotor
        FemmProblem=rotategroups_mfemm(FemmProblem,2,th_m);
        %mi_selectgroup(2), mi_moverotate(0,0,th_m);
        % redraw the airgap arc
        FemmProblem=draw_airgap_arc_with_meshX(FemmProblem,geo,th_m,geo.fem,boundnameAPmove);
        
        %mi_saveas([pathname,'\mot_temp.fem']);
        plotfemmproblem(FemmProblem);
        filename = 'mot0.fem';
        writefemmfile(filename, FemmProblem);
        
        %filename = fmesher(filename);
        %ansfile = fsolver(filename);
        %myfpproc = fpproc();
        %myfpproc.opendocument(ansfile);
        
        %mi_analyze(1);
        
        %mi_loadsolution;
        %         keyboard
        %post_proc;
        %mo_close, mi_close
        
        SOL = [SOL; sol];
    end
    % ripple_pu = abs(std(SOL(:,6))/mean(SOL(:,6)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       FINE VECCHIO CICLO FOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%% VERIFICA CHE I DUE CICLI FOR E PARFOR DIANO LO STESSO RISULTATO %%
% ERROR = SOL-SOL_PARFOR_temp
