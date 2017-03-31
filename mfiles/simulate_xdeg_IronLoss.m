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

function [SOL] = simulate_xdeg_IronLoss(geo,io,Br,gamma_in,eval_type,mat,per)

% pathname = geo.pathname;
pathname = cd;

th0 = geo.th0;
p   = geo.p;
q = geo.q;
r  = geo.r;
gap = geo.g;
ns  = geo.ns;
pc  = 360/(ns*p)/2;     % span of half slot [mach degree]
ps  = geo.ps;
Qs = geo.Qs;
l = geo.l;

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
        xdeg = max(60,20*Qs/q);     % rotation excursion is 60 ele degree or one third of pole pitch
        nsim = max(geo.nsim_singt,6);
        gamma = gamma_in;
        randFactor = 0;
    case 'singm'
        xdeg = max(60,20*Qs/q);
        nsim = max(geo.nsim_singt,6);
        gamma = gamma_in;
        randFactor = 0;
end

geo.nsim_singt = nsim; geo.delta_sim_singt = xdeg;

%% simulation angle
gradi_da_sim=180/p*ps;

id = io * cos(gamma * pi/180);
iq = io * sin(gamma * pi/180);

Hc = Br/(4e-7*pi);

SOL = [];
FluxDens = [];

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

openfemm
%main_minimize
opendocument([pathname,'\mot0.fem']);

for jj = 1:nsim
    
    th_m = (th(jj) - th0)/p;
%     openfemm
%     opendocument([pathname,'\mot0.fem']);
    
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
    if strcmp(geo.RotType,'SPM')
        mi_modifymaterial('Bonded-Magnet',3,Hc);
    else
        if length(Hc)==1
            Hc_vect = Hc*ones(1,length(geo.BLKLABELS.rotore.BarName));
        else
            Hc_vect=[Hc Hc];
        end
        for ii = 1:length(Hc_vect)
            mi_modifymaterial(['Bonded-Magnet' num2str(ii)],3,Hc_vect(ii));
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
    post_proc_IronLoss;
%     mo_close, mi_close
%     closefemm
    
    SOL = [SOL; sol];
end
mo_close, mi_close
closefemm

% Effective value of flux and current, simulation are done with one turns
% in slot and consequently, current in fem simulation is increase by the number of conductors in slot Nbob....
SOL(:,2)=SOL(:,2)/geo.Nbob;
SOL(:,3)=SOL(:,3)/geo.Nbob;
SOL(:,4)=SOL(:,4)*geo.Nbob;
SOL(:,5)=SOL(:,5)*geo.Nbob;

if io == 0 && strcmp(geo.BLKLABELSmaterials(6),'Air')
    temp = zeros(nsim-4,1);
    Pfes_h = 0;
    Pfes_c = 0;
    Pfer_h = 0;
    Pfer_c = 0;
    SOL(:,7) = [Pfes_h;Pfes_c;Pfer_h;Pfer_c;temp];
else
    %% calculate losses
    FluxDens = [FluxDens; fluxdens];
    
    % Iron Factor (stator)
    Kh_s = mat.Stator.kh;
    alpha_s = mat.Stator.alpha;
    beta_s = mat.Stator.beta;
    Ke_s = mat.Stator.ke;
    MassDensity_s = mat.Stator.kgm3;                     % Mass density [kg/m^3]
    
    % Iron Factor (rotor)
    Kh_r = mat.Rotor.kh;
    alpha_r = mat.Rotor.alpha;
    beta_r = mat.Rotor.beta;
    Ke_r = mat.Rotor.ke;
    MassDensity_r = mat.Rotor.kgm3;                     % Mass density [kg/m^3]
    
    EvalSpeed = per.EvalSpeed;
    
    Fw = p*EvalSpeed/60;                           % Working Frequency [Hz],
    T = 1/Fw;                                      % rotation period
    Fs = nsim/T;                                   % sampling frequency
    Ts = 1/Fs;                                     % sampling period
    f = Fs*(0:(nsim/2))/nsim;                      % Real frequency
    t = (0:nsim-1)*Ts;                             % Real time vector
    
    %% extrat data
    posIron = FluxDens(:,1);            % Element postion
    areaIron = FluxDens(:,2);            % Iron area  [mm^2]
    groNoIron = FluxDens(:,3);          % Element group (12 Stator iron - 22 Rotor iron)
    
    % Recover flux density info
    fIron = FluxDens(:,4:end);
    fIron = fIron.';
    
    %% stator iron block flux density
    [line, row] = find(groNoIron == 12);    % select elements belongs to stator iron
    
    for ii = 1:size(line,1)
        fIronSta(:,ii) = fIron (:,line(ii));
        areaIronSta(ii,:) = areaIron(line(ii),:);
        posIronSta(ii,:) = posIron(line(ii),:);
    end
    
    fIronSta(find(fIronSta(:,2)==0),:)=[];
    posAngle = angle(posIronSta);       % angle from xoy axe
    
    B_Sta_x =real(fIronSta);            % Bx
    B_Sta_y =imag(fIronSta);            % By
    
    %% XY to RT for stator iron  A_180 = [A_60+(-C_60)+(-B_60)]
    for jj = 1:size(posAngle)
        for ii = 1:size(B_Sta_x,1)
            if ii <= size(B_Sta_x,1)/3
                B_Sta_r(ii,jj) = B_Sta_x(ii,jj)*cos(posAngle(jj)) + B_Sta_y(ii,jj)*sin(posAngle(jj));
                B_Sta_t(ii,jj) = B_Sta_y(ii,jj)*cos(posAngle(jj)) - B_Sta_x(ii,jj)*sin(posAngle(jj));
            elseif ii > size(B_Sta_x,1)/3 && ii <= size(B_Sta_x,1)*2/3
                B_Sta_r(ii,jj) = B_Sta_x(ii,jj)*cos(posAngle(jj)+60/p/180*pi) + B_Sta_y(ii,jj)*sin(posAngle(jj)+60/p/180*pi);
                B_Sta_t(ii,jj) = B_Sta_y(ii,jj)*cos(posAngle(jj)+60/p/180*pi) - B_Sta_x(ii,jj)*sin(posAngle(jj)+60/p/180*pi);
                B_Sta_t(ii,jj) = - B_Sta_t(ii,jj);
                B_Sta_r(ii,jj) = - B_Sta_r(ii,jj);
            else
                B_Sta_r(ii,jj) = B_Sta_x(ii,jj)*cos(posAngle(jj)+2*60/p/180*pi) + B_Sta_y(ii,jj)*sin(posAngle(jj)+2*60/p/180*pi);
                B_Sta_t(ii,jj) = B_Sta_y(ii,jj)*cos(posAngle(jj)+2*60/p/180*pi) - B_Sta_x(ii,jj)*sin(posAngle(jj)+2*60/p/180*pi);
                B_Sta_t(ii,jj) = - B_Sta_t(ii,jj);
                B_Sta_r(ii,jj) = - B_Sta_r(ii,jj);
            end
        end
    end
    
    %% transfer phase B & C data to phase A, A_180 = [A_60+(-C_60)+(-B_60)]
    temp = B_Sta_t(6+1:2*6,:);
    B_Sta_t(6+1:2*6,:) = B_Sta_t(2*6+1:end,:);
    B_Sta_t(2*6+1:end,:) = temp;
    %% 180 ele degree to 360 ele degree
    B_Sta_t = [B_Sta_t; -B_Sta_t];
    temp = B_Sta_r(6+1:2*6,:);
    B_Sta_r(6+1:2*6,:) = B_Sta_r(2*6+1:end,:);
    B_Sta_r(2*6+1:end,:) = temp;
    %% 180 ele degree to 360 ele degree
    B_Sta_r = [B_Sta_r; -B_Sta_r];
    
    for jj = 1:size(line)
        [B_FFT_hoSta_r(:,jj),B_FFT_magSta_r(:,jj)] = FFTAnalysis(B_Sta_r(:,jj),Fw);
        [B_FFT_hoSta_t(:,jj),B_FFT_magSta_t(:,jj)] = FFTAnalysis(B_Sta_t(:,jj),Fw);
    end
    
    B_FFT_magSta = sqrt(B_FFT_magSta_r.^2+B_FFT_magSta_t.^2);
    
    PhystStapu = 3/2*Kh_s* (B_FFT_hoSta_r(2,:).^alpha_s).*B_FFT_magSta(2,:).^beta_s;              % [W/kg]
    PeddyStapu = Ke_s* (B_FFT_hoSta_r(1:17,:).*B_FFT_magSta(1:17,:)).^2;                    % [W/kg]
    
    for ii =1:length(PhystStapu)
        PhystSta(:,ii) = PhystStapu(:,ii) * areaIronSta(ii) * l * MassDensity_s*1e-9;                           % [W]
        PeddySta(:,ii) = PeddyStapu(:,ii) * areaIronSta(ii) * l * MassDensity_s*1e-9;                           % [W]
    end
    
    %% Stator iron loss
    Pfes_h = sum(sum(PhystSta)) * 2*p/ps *3;        % whole stator hysteresis loss
    Pfes_c = sum(sum(PeddySta)) * 2*p/ps *3;        % whole stator eddy current loss
    
    %% extract rotor iron block flux density
    [line, row] = find(groNoIron == 22);    % select elements belongs to stator iron
    
    for ii = 1:size(line,1)
        fIronRot(:,ii) = fIron (:,line(ii));
        areaIronRot(ii,:) = areaIron(line(ii),:);
        posIronRot(ii,:) = posIron(line(ii),:);
    end
    
    fIronRot(find(fIronRot(:,2)==0),:)=[];
    
    B_Rot_d = real(fIronRot);
    B_Rot_q = imag(fIronRot);
    
    %% recover data for 360 ele degree rotation
    for jj = 1:size(line)
        B_Rot_d_360(:,jj) = repmat(B_Rot_d(:,jj),p*3,1);
        B_Rot_q_360(:,jj) = repmat(B_Rot_q(:,jj),p*3,1);
    end
    
    for jj = 1:size(line)
        [B_FFT_hoRot_d(:,jj),B_FFT_magRot_d(:,jj)] = FFTAnalysis(B_Rot_d_360(:,jj),Fw);
        [B_FFT_hoRot_q(:,jj),B_FFT_magRot_q(:,jj)] = FFTAnalysis(B_Rot_q_360(:,jj),Fw);
    end
    
    B_FFT_magRot = sqrt(B_FFT_magRot_d.^2+B_FFT_magRot_q.^2);
    
    PhystRotpu = 3/2 * Kh_r*((B_FFT_hoRot_d(1,:)+Fw).^alpha_r).*B_FFT_magRot(1,:).^beta_r;
    PeddyRotpu = Ke_r* ((B_FFT_hoRot_d(1:19,:)+Fw).*B_FFT_magRot(1:19,:)).^2;
    
    for ii =1:length(PeddyRotpu)
        PhystRot(:,ii) = PhystRotpu(:,ii) * areaIronRot(ii) * l * MassDensity_r*1e-9;                           % [W]
        PeddyRot(:,ii) = PeddyRotpu(:,ii) * areaIronRot(ii) * l * MassDensity_r*1e-9;                           % [W]
    end
    
    Pfer_h = sum(sum(PhystRot))* 2*p/ps;        % whole rotor hysteresis loss
    Pfer_c = sum(sum(PeddyRot))* 2*p/ps;        % whole rotor eddy current loss
    
    %% save losses data
    temp = zeros(nsim-4,1);
    SOL(:,7) = [Pfes_h;Pfes_c;Pfer_h;Pfer_c;temp];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       FINE VECCHIO CICLO FOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% VERIFICA CHE I DUE CICLI FOR E PARFOR DIANO LO STESSO RISULTATO %%
% ERROR = SOL-SOL_PARFOR_temp
