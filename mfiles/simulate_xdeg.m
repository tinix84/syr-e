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

% - Output structure: SOL

function [SOL] = simulate_xdeg(geo,per,mat,eval_type,pathname,filename,LossEvaluationCheck)

% NB: pathname with final slash
% eval type determines number of slimulated positions

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
        if (geo.periodicity==4)
            xdeg = min(geo.delta_sim_singt,120);
        end
        if (geo.periodicity==5)
            xdeg = min(geo.delta_sim_singt,60);
        end
        nsim = round(geo.nsim_singt*xdeg/geo.delta_sim_singt);
        gamma = per.gamma;
        randFactor = 0;
    case 'singm'
        if (geo.periodicity == 4)
            xdeg = min(geo.delta_sim_singt,120);
        end
        if (geo.periodicity == 5)
            xdeg = min(geo.delta_sim_singt,60);
        end
        nsim = round(geo.nsim_singt*xdeg/geo.delta_sim_singt);
        gamma = per.gamma;
        randFactor = 0;
end

th0 = geo.th0;
p   = geo.p;
% r  = geo.r;
% gap = geo.g;
ns  = geo.ns;
pc  = 360/(ns*p)/2;
ps  = geo.ps;
n3phase = geo.n3phase; %AS number of 3-phase circuits

% rotor positions
if strcmp(eval_type,'MO_OA')||strcmp(eval_type,'MO_GA')
    % during optimization, random position offset
    sim_step=xdeg/(nsim+0.5);
    offset=1*sim_step*rand;
else
    % during re-evaluation, regular position steps
    sim_step=xdeg/(nsim);
    offset=0;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% ciclo for %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open and draw motor once, rotate and simulate nsim positions by Chao Lu 17/01/2017
openfemm
opendocument([pathname,filename]);

SOL.th = zeros(1,nsim); % electrical angle in degree
SOL.id = zeros(1,nsim); % d-axis current
SOL.iq = zeros(1,nsim); % q-axis current
SOL.fd = zeros(1,nsim); % d-axis flux linkage
SOL.fq = zeros(1,nsim); % q-axis flux linkage
SOL.T  = zeros(1,nsim); % Torque from block integral (suggested by David Meeker)
SOL.F  = zeros(1,nsim); % radial and tangential force of the simulated portion

phase_name = cell(n3phase*3,1);
phase_name_neg = cell(n3phase*3,1);
for jj = 1:nsim
    
    th_m = (th(jj) - th0)/p;
    
    % assign the phase current values to the FEMM circuits
    for ik=0:(n3phase-1)
        
        if geo.avv_flag((3*ik)+1)==1 && geo.avv_flag((3*ik)+2)==1 && geo.avv_flag((3*ik)+3)==1
            i123 = dq2abc(id,iq,th(jj)*pi/180,n3phase,ik);
            i_tmp((3*ik)+1,jj) = i123(1);
            i_tmp((3*ik)+2,jj) = i123(2);
            i_tmp((3*ik)+3,jj) = i123(3);
        else
            if geo.avv_flag((3*ik)+1)==0 && geo.avv_flag((3*ik)+2)==0 && geo.avv_flag((3*ik)+3)==0
                i_tmp((3*ik)+1,jj) = 0;
                i_tmp((3*ik)+2,jj) = 0;
                i_tmp((3*ik)+3,jj) = 0;
            end
        end
        
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
        mi_modifymaterial(mat.LayerMag.MatName,3,mat.LayerMag.Hc);
    else
        if length(mat.LayerMag.Hc)==1
            Hc_vect = mat.LayerMag.Hc*ones(1,length(geo.BLKLABELS.rotore.BarName));
        else
            %             Hc_vect=[Hc Hc];  // DUBBIO -- a cosa serve? 2018 07 26
        end
        for ii = 1:length(Hc_vect)
            mi_modifymaterial([mat.LayerMag.MatName '_' num2str(ii)],3,Hc_vect(ii));
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
    
    %%     post_proc  WAS HERE
    
    % load phase flux linkages
    for ii=0:(n3phase-1) %AS
        temp_out = mo_getcircuitproperties(phase_name{3*ii+1});
        temp_out = temp_out - mo_getcircuitproperties(phase_name_neg{3*ii+1});
        f(3*ii+1) = temp_out(3) * 2 * p/ps; %ps number of poles in FEMM
        temp_out = mo_getcircuitproperties(phase_name{3*ii+2});
        temp_out = temp_out - mo_getcircuitproperties(phase_name_neg{3*ii+2});
        f(3*ii+2) = temp_out(3) * 2 * p/ps;
        temp_out = mo_getcircuitproperties(phase_name{3*ii+3});
        temp_out = temp_out - mo_getcircuitproperties(phase_name_neg{3*ii+3});
        f(3*ii+3) = temp_out(3) * 2 * p/ps;
    end
    
    for ik=0:(n3phase-1) %AS
        fdq = abc2dq(f(3*ik+1),f(3*ik+2),f(3*ik+3),th(jj)*pi/180,n3phase,ik);
        fd_temp(ik+1,jj)=fdq(1);
        fq_temp(ik+1,jj)=fdq(2);
    end
    
    fd=mean(fd_temp(:,jj));
    fq=mean(fq_temp(:,jj));
    
    % block evaluation of torque and force
    for ii=1:length(geo.BLKLABELS.rotore.xy(:,1))
        xB=geo.BLKLABELS.rotore.xy(ii,1);
        yB=geo.BLKLABELS.rotore.xy(ii,2);
        [xB,yB]=rot_point(xB,yB,th_m*pi/180);
        mo_selectblock(xB,yB);
    end
    
    Tblock=mo_blockintegral(22)*2*p/ps;
    mo_clearblock;
    
    % Calcolo del Volume dei magneti VolPM - rev.Gallo 14/03/2018
    if jj==1 %viene calcolato solo alla prima simulazione (rotore in posizione di partenza)
        flagPM=0;
        for ii=1:length(geo.BLKLABELS.rotore.xy(:,1))
            if geo.BLKLABELS.rotore.xy(ii,3)== 6 %lettura codice del materiale per individuare regioni di PM presenti nella struttura di rotore
                xA=geo.BLKLABELS.rotore.xy(ii,1);
                yA=geo.BLKLABELS.rotore.xy(ii,2);
                [xA,yA]=rot_point(xA,yA,th_m*pi/180); %rotazione delle coordinate di random position offset
                mo_selectblock(xA,yA);
                flagPM=1;
            end
        end
        if flagPM %flag per capire se area del magnete rettangolare è presente o no
            VolPM=(2*geo.p*mo_blockintegral(10))/geo.ps; %Calcolo Volume magnete totale nel rotore [m3]
        else
            VolPM=0;
        end
        
        mo_clearblock();
    end
    
    
    
    %%     post_proc  WAS HERE - end
    
    SOL.th(jj) = th(jj);
    SOL.id(jj) = id;
    SOL.iq(jj) = iq;
    SOL.fd(jj) = fd;  %AS
    SOL.fq(jj) = fq;  %AS
    SOL.T(jj)  = Tblock;
    SOL.VolPM(jj)= VolPM;
    
end

if(LossEvaluationCheck)
    % if jj == 1
    if jj == nsim
        EleNo = mo_numelements;               % Number of mesh elements
        pos = zeros(EleNo,1);                 % Matrix that will hold the mesh elements centroid coordinates as complex number
        area = zeros(EleNo,1);                % Matrix that will hold the mesh elements area
        groNo = zeros(EleNo,1);               % Matrix that will hold the mesh elements group number
        for i = 1:EleNo
            elm = mo_getelement(i);
            pos(i) = elm(4)+j*elm(5);
            area(i) = elm(6);
            groNo(i) = elm(7);
        end
        
        %% find stator iron elements
        [line_Sta, row] = find(groNo==12);      % Stator Iron
        pos_Sta = pos(line_Sta);
        area_Sta = area(line_Sta);
        gro_Sta = groNo(line_Sta);
        %% select 1/3 stator iron elements
        Sta = find(angle(pos_Sta) < ((geo.Qs/3*2-1)*pc*pi/180));
        pos_Sta = pos_Sta(Sta);
        area_Sta = area_Sta(Sta);
        gro_Sta = gro_Sta(Sta);
        
        %% find rotor iron element
        [line_Rot, row] = find(groNo==22);
        pos_Rot = pos(line_Rot);
        area_Rot = area(line_Rot);
        gro_Rot = groNo(line_Rot);
        
        %% Matrix to store data
        areaIron = [area_Sta; area_Rot];
        posIron = [pos_Sta; pos_Rot];
        groNoIron = [gro_Sta; gro_Rot];
        %% delete redundant data, just keep iron element data
        IronNo = size(groNoIron,1);
        fIron =zeros(IronNo,6);            % Matrix that save flux density of iron element
    end
    
    for kk = 1:1:5
        if jj == 1+(kk-1)*round(nsim/5)               % 5 positions are simulated for losses
            RotPos = exp(j*th_m*pi/180);              % It is a parameter that considers the rotor position
            for i = 1:IronNo
                if (groNoIron(i) == 22)                       % rotor iron elements
                    Pos_Rot = posIron(i)*RotPos;              % get element position
                    %% axis transform
                    fIron(i,kk) = (mo_getb(real(Pos_Rot),imag(Pos_Rot))*[1;j]);
                    Bd = real(fIron(i,kk)) * cosd(th_m) + imag(fIron(i,kk)) * sind(th_m);
                    Bq = imag(fIron(i,kk)) * cosd(th_m) - real(fIron(i,kk)) * sind(th_m);
                    fIron(i,kk) = [Bd,Bq] * [1;j];
                elseif (groNoIron(i) == 12)                   % stator iron elements
                    if ceil(th_m) < 60/p                      % for stator, just 60 ele degree rotation data is needed
                        %%          first one third of stator
                        Pos_Sta = posIron(i);
                        fIron(i,kk) = (mo_getb(real(Pos_Sta),imag(Pos_Sta))*[1;j]);
                        %%          second one third stator iron area
                        Pos_Sta1 = Pos_Sta*exp(j*(geo.Qs/3*2*pc)*pi/180);
                        fIron(i,kk+6) = (mo_getb(real(Pos_Sta1),imag(Pos_Sta1))*[1;j]);
                        %%          last one third stator iron area
                        Pos_Sta2 = Pos_Sta1*exp(j*(geo.Qs/3*2*pc)*pi/180);
                        fIron(i,kk+6*2) = (mo_getb(real(Pos_Sta2),imag(Pos_Sta2))*[1;j]);
                    end
                end
            end
        end
    end
    
    if jj == nsim                                 % Last postion data collection
        RotPos = exp(j*th_m*pi/180);              % It is a parameter that considers the rotor position
        for i = 1:IronNo
            if (groNoIron(i) == 22)                       % rotor iron elements
                Pos_Rot = posIron(i)*RotPos;              % get element position
                %% axis transform
                fIron(i,6) = (mo_getb(real(Pos_Rot),imag(Pos_Rot))*[1;j]);
                Bd = real(fIron(i,6)) * cosd(th_m) + imag(fIron(i,6)) * sind(th_m);
                Bq = imag(fIron(i,6)) * cosd(th_m) - real(fIron(i,6)) * sind(th_m);
                fIron(i,6) = [Bd,Bq] * [1;j];
            elseif (groNoIron(i) == 12)                   % stator iron elements
                if th_m < 60/p                      % for stator, just 60 ele degree rotation data is needed
                    %%          first one third of stator
                    Pos_Sta = posIron(i);
                    fIron(i,6) = (mo_getb(real(Pos_Sta),imag(Pos_Sta))*[1;j]);
                    %%          second one third stator iron area
                    Pos_Sta1 = Pos_Sta*exp(j*(geo.Qs/3*2*pc)*pi/180);
                    fIron(i,6+6) = (mo_getb(real(Pos_Sta1),imag(Pos_Sta1))*[1;j]);
                    %%          last one third stator iron area
                    Pos_Sta2 = Pos_Sta1*exp(j*(geo.Qs/3*2*pc)*pi/180);
                    fIron(i,6+6*2) = (mo_getb(real(Pos_Sta2),imag(Pos_Sta2))*[1;j]);
                end
            end
        end
        fluxdens = [posIron areaIron groNoIron fIron];      % Saving all the info necessary to calculate Iron Losses
    end
    
    
end


if(LossEvaluationCheck)
    %         if io == 0 && strcmp(geo.BLKLABELSmaterials(6),'Air')
    %             %temp = zeros(nsim-4,1);
    %             Pfes_h = 0;
    %             Pfes_c = 0;
    %             Pfer_h = 0;
    %             Pfer_c = 0;
    %             %SOL(:,7) = [Pfes_h;Pfes_c;Pfer_h;Pfer_c;temp];
    %             SOL.Pfes_h = Pfes_h;
    %             SOL.Pfes_c = Pfes_c;
    %             SOL.Pfer_h = Pfer_h;
    %             SOL.Pfer_c = Pfer_c;
    %         else
    %% calculate losses
    FluxDens = fluxdens;
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
        fIronSta(:,ii) = fIron(:,line(ii));
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
    %     temp = B_Sta_t(6+1:2*6,:);
    %     B_Sta_t(6+1:2*6,:) = B_Sta_t(2*6+1:end,:);
    %     B_Sta_t(2*6+1:end,:) = temp;
    
    temp1 = size(B_Sta_t,1)/3;
    
    temp = B_Sta_t(temp1+1:2*temp1,:);
    B_Sta_t(temp1+1:2*temp1,:) = B_Sta_t(2*temp1+1:end,:);
    B_Sta_t(2*temp1+1:end,:) = temp;
    %% 180 ele degree to 360 ele degree
    B_Sta_t = [B_Sta_t; -B_Sta_t];
    
    %     temp = B_Sta_r(6+1:2*6,:);
    %     B_Sta_r(6+1:2*6,:) = B_Sta_r(2*6+1:end,:);
    %     B_Sta_r(2*6+1:end,:) = temp;
    temp = B_Sta_r(temp1+1:2*temp1,:);
    B_Sta_r(temp1+1:2*temp1,:) = B_Sta_r(2*temp1+1:end,:);
    B_Sta_r(2*temp1+1:end,:) = temp;
    %% 180 ele degree to 360 ele degree
    B_Sta_r = [B_Sta_r; -B_Sta_r];
    
    for jj = 1:size(line)
        [B_FFT_hoSta_r(:,jj),B_FFT_magSta_r(:,jj)] = FFTAnalysis(B_Sta_r(:,jj),Fw);
        [B_FFT_hoSta_t(:,jj),B_FFT_magSta_t(:,jj)] = FFTAnalysis(B_Sta_t(:,jj),Fw);
    end
    
    B_FFT_magSta = sqrt(B_FFT_magSta_r.^2+B_FFT_magSta_t.^2);
    
    PhystStapu = 3/2*Kh_s* (B_FFT_hoSta_r(2,:).^alpha_s).*B_FFT_magSta(2,:).^beta_s; % [W/kg]
    PeddyStapu = Ke_s* (B_FFT_hoSta_r(1:17,:).*B_FFT_magSta(1:17,:)).^2;             % [W/kg]
    
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
    %temp = zeros(nsim-4,1);
    %SOL(:,7) = [Pfes_h;Pfes_c;Pfer_h;Pfer_c;temp];
    SOL.Pfes_h = Pfes_h;
    SOL.Pfes_c = Pfes_c;
    SOL.Pfer_h = Pfer_h;
    SOL.Pfer_c = Pfer_c;
    %         end
    
end

mo_close, mi_close
closefemm

% Apply number of turns (simulation done with one turn per coil)
SOL.id=SOL.id/geo.Nbob/n3phase;
SOL.iq=SOL.iq/geo.Nbob/n3phase;
SOL.fd=SOL.fd*geo.Nbob*n3phase;
SOL.fq=SOL.fq*geo.Nbob*n3phase;

