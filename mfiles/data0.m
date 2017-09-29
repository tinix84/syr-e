function [bounds, objs, geo, per, mat] = data0(dataIn)

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data0.m:
% manual input data used as default by the graphic input GUI
% used also with the GUI

% per: performance
% geo: geometry
% bounds: bounds of optimization inputs RQ

if (nargin<1)
    
    %% MANUAL SETTINGS
    % performance target
    per.Loss = 1000;                % Joule loss [W]
    per.tempcu = 130;               % Target Copper Temperature [C]
    per.overload = 2;               % current overload factor used for optimization (1 means Joule loss = per.Loss)
    per.BrPP = 0.4;                 % Br used for postprocessing [T]
    per.temphous = 70;              % Housing Temperature [C]
%     per.tempcuest = per.tempcu;     % Estimated Copper Temperature [C]
    % torque and ripple penalization
    % a penalizing coefficient will be applied during the optimization to all
    % the machines with a lower torque or a higher torque ripple
    per.min_exp_torque = 10;   % minimum expected torque [Nm]
    per.max_exp_ripple = 8;    % maximum expected torque ripple in % of T
    per.max_Cu_mass = 10;       % maximum expected copper mass [kg]
    per.io = 0;
    
    %% geo: geometry
    geo.BLKLABELSmaterials={
        'Air';
        'Air';
        'Copper';
        '10JNEX900';
        '10JNEX900';
        'NdFeB 37 MGOe';
        '10JNEX900';
        'Copper'};
    
    % create dataSet only for material definition
    dataIn.SlotMaterial = geo.BLKLABELSmaterials{3};
    dataIn.StatorMaterial = geo.BLKLABELSmaterials{4};
    dataIn.RotorMaterial = geo.BLKLABELSmaterials{5};
    dataIn.FluxBarrierMaterial = geo.BLKLABELSmaterials{6};
    dataIn.ShaftMaterial = geo.BLKLABELSmaterials{7};
    dataIn.RotorCondMaterial = geo.BLKLABELSmaterials{8};
    
    geo.pont0 = 0.4;    % thickness of the structural bridges at the airgap
    % AND minimum mechanical tolerance [mm]
    
    % Type of barriers
    % 'Circular': Circular barriers
    % 'ISeg'    : I-shaped exterior barrier plus Seg(mented) or U-shaped barriers
    % 'Seg'     : All Seg(mented) barriers
    % 'Fluid'   : barriers shaped according to fluid flow lines
    % 'SPM'     : Surface mounted permanent magnet motor
    geo.RotType = 'Circular';
    
    geo.RemoveTMPfile = 'ON';      % 'ON' of 'OFF' for remuving the motor folders in tmp
    
    geo.RaccBarrier='=OFF';  % Seg geometry, it attivate circular connection to the side of the barrier if 'ON'
    % Seg geometry, barrier are drawn with no circular junction og if 'OFF'
    geo.DTrasl=0;
    
    geo.p  = 2;         % pole pairs
    geo.R  = 100;      % stator outer radius [mm]
    geo.r  = 0.65*geo.R;%0.6*geo.R-0.275;  % airgap radius [mm]
    geo.Ar = 25;     % shaft radius [mm]
    geo.g  = 0.5;     % airgap [mm]
    geo.l  = 215;       % stack length [mm]
    
    geo.BarFillFac=1;   % barrier filling factor
    % stator
    geo.q   = 1;        % stator slots per pole per phase
    %     geo.slot_layer_pos='over_under';  %     stator slot layer position, two
    %     solution are possible, 'over_under', and 'side_by_side'
    geo.slot_layer_pos='side_by_side';
    % tooth length and width
    geo.b  = 0.5;           % Bgap/Bfe,yoke  (back-iron p.u. size)
    geo.bt = 1;           % Bgap/Bfe,tooth (tooth p.u. size)
    
    geo.lt = geo.R * (1-geo.r/geo.R*(1+geo.b/geo.p));    % tooth length
    geo.wt = 5;  %pi/(3*geo.q*geo.p)*geo.r*geo.bt;       % tooth width
    
    geo.acs  = 0.25;    % stator slot opening [p.u.]
    geo.ttd = 1;        % tooth tang depth [mm]
    geo.tta = 10;       % tooth tang angle (mech degree)
    geo.SFR = 1;        % fillet at the back corner of the slot [mm]
    
    % rotor
    if strcmp(geo.RotType,'SPM')
        geo.naly = 1;
    else
        geo.nlay  = 3;      % number of layers
    end
    
    geo.racc_pont = 1 * geo.pont0;                  % radius of the fillet at the sides of inner bridges (if any) [mm]
    geo.ang_pont0 = geo.pont0 / geo.r * 180/pi;     % span of the arc corresponding to pont0 at the airgap radius [deg]
    geo.hfe_min = 2*geo.pont0;                      % minimum tickness of each steel flux guide
    % winding description
    geo.kcu = 0.44;                                % slot filling factor (net copper/slot area)
    %    geo.avv   = [1 1 -3 -3 2 2                      % winding scheme
    %                 1 1 -3 -3 2 2];
    geo.avv=[1 2 3 ;
        -3 -1 -2];
    
    geo.kracc = 9/9;                                % pitch shortening factor (for end connections length estimation)
    geo.Ns    = 87;                                % turns in series per phase (entire motor, one way, all poles in series)
    geo.ns = geo.q*6;                               % number of slot per pole pair
    geo.Nbob  = geo.Ns/geo.p/(geo.ns/6)/size(geo.avv,1);  % conductors in slot per label
    geo.nmax = 4000;                                % overspeed [rpm]
    geo.Qs = 3;                                     % number of slots in the simulation
    
    % permanent magnets (default: Br = 0)
%     geo.Br = 0.0;
%     geo.Hc = 1/(4e-7*pi)*geo.Br;
%     geo.rhoPM = 7600;                               % Permanent Magnet mass density [Kg/m3] (for NdFeB magnets)
%     geo.Br_commercial=1.25;                         % [T] remanence of commercial NdFeB permanent magnets
    
    geo.lm = 5;      % the thickness of permanent magnet
    geo.phi = 5/6*180;       % the angle range of permanent magnet
    % direction of magnetization in PMs of SPM with multiple segments of PM
    geo.PMdir = 'p';    % parallel direction
    geo.PMdir = 'r';    % radial direction
    
    % Settings that vary if MOOA or not (MOOA means during the
    % optimization)
    % Mesh ratio (all_motor/air-gap)
    geo.K_mesh_MOOA = 10;    % optimization
    geo.K_mesh = 2;         % post-processing and manual design
    
    % number of simulated rotor positions
    geo.nsim_MOOA = 5;          % simulated positions
    geo.randFactor = 0;         % Noise factor for position number reduction
    geo.delta_sim_MOOA = 360/(geo.q*6);    % rotor position span [elt degrees]
    % evalx means the re-evaluation stage, with a finer resolution
    geo.nsim_singt = 30;        % simulated positions
    geo.delta_sim_singt = 60;   % rotor position span [elt degrees]
    
    % rotor geometry (subject to optimization: default set to zero)
    geo.x0=[0 0]; % geo.alpha = 0;
    geo.x0 = geo.r/cos(pi/2/geo.p);
    geo.dalpha_pu = [0.45 0.22*ones(1,geo.nlay-1)];
    geo.hc_pu = 0.5*ones(1,geo.nlay);
    if strcmp(geo.RotType,'SPM')
        geo.dx = 1;
    else
        geo.dx = zeros(1,geo.nlay);
    end
    
    geo.radial_ribs_eval = 0;   % if 1 the radial ribs are manually designed, if 0, the radial ribs are automatic designed
    geo.pont = [0 0 0];         % radial ribs
    
    %% bounds: limits of the search space
    % dalpha1 [p.u.]
    bounds_dalpha_1 = [22.5 45]/90;   % first angle [deg]
    % dalpha(2:nlay) [p.u.]
    bounds_dalpha = ones(geo.nlay-1,1) * [0.5/geo.nlay 0.5];% other angles [p.u.]
    if strcmp(geo.RotType, 'SPM')
        % thickness of PM
        bounds_hc = geo.lm * [0.7 1.5];
        geo.hc_pu = mean(bounds_hc);
    else
        % barrier ticknesses [p.u.]
        bounds_hc = ones(geo.nlay,1) * [0.2 1];
    end
    % barriers radial offset [p.u.]
    bounds_dx = ones(geo.nlay,1) * [-0.75 0.75];
    % remanence of the PMs
    bounds_Br = ones(geo.nlay,1) * [0.2 0.8];
    % airgap
    bounds_g = geo.g * [0.7  1.5];
    % Rotor radius
    bounds_xr = geo.r * [0.8  1.2];
    % Tooth width [mm]
    bounds_wt = geo.wt * [0.75  1.25];
    % Tooth length [mm]
    bounds_lt = geo.lt * [0.8 1.2];
    %  stator slot opening [p.u.]
    bounds_acs = [0.15 0.5];
    % tooth tang depth [mm]
    bounds_ttd = geo.ttd * [0.8 1.2];
    % phase angle of the current vector
    bounds_gamma = [40 75];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% RQ and bounds: define here the variables to be optimized
    %% RQ is the vector of the optimization inputs (variable size n)
    %% bounds is the n x 2 vector containig the boundaries of each input
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yes_vector = ones(geo.nlay,1);
    no_vector = zeros(geo.nlay,1);
    if (strcmp(geo.RotType,'Fluid') || strcmp(geo.RotType,'Seg') || strcmp(geo.RotType,'Circular'))
        flag_dx = yes_vector;
    else
        flag_dx = no_vector;
    end
    
    % bounds
    % put 1 or 0 to add/remove variables from RQ and bounds
    % example:
    if geo.nlay == 1
        bounds = [
            bounds_dalpha_1 1
            bounds_dalpha   yes_vector(1)
            bounds_hc       1
            bounds_dx       0
            bounds_Br       0
            bounds_g        0
            bounds_xr       0
            bounds_wt       0
            bounds_lt       0
            bounds_acs      0
            bounds_ttd      0
            bounds_gamma    0];
    else
        bounds = [
            bounds_dalpha_1 1
            bounds_dalpha   yes_vector(1:end-1)
            bounds_hc       yes_vector
            bounds_dx       yes_vector*1
            bounds_Br       yes_vector*0
            bounds_g        0
            bounds_xr       0
            bounds_wt       0
            bounds_lt       0
            bounds_acs      0
            bounds_ttd      0
            bounds_gamma    0];
    end
    
    % eliminate unnecessary rows
    filt_bounds = (bounds(:,3)==1);
    bounds = bounds(filt_bounds,1:2);
    bounds = round(bounds*100)/100;
    %     % eliminate the unnecessary RQnames
    %     RQnames = RQnames(filt_bounds);
    %     geo.RQnames = RQnames;
    
    %% OBJECTIVES
    objs = [per.min_exp_torque  1      % minimum expected torque
            per.max_exp_ripple  1      % maximum expected torque ripple
            per.max_Cu_mass     0];    % maximum expected copper mass
    filt_objs = (objs(:,2)==1);
    objs = objs(objs(:,2)==1,:);
    per.objs=objs;

    
else
    
    %% READ INPUTS FROM THE GUI
    
    % main performance target
    per.Loss = dataIn.AdmiJouleLosses;             % admitted Joule loss [W]
    per.tempcu = dataIn.TargetCopperTemp;           % Target Copper Temperature [C]
    %     per.Vdc = dataIn.DCVoltage;              % dc-link voltage [V]
    per.overload = dataIn.CurrOverLoad;           % current overload factor used for optimization (1 means Joule loss = per.Loss)
    per.BrPP = dataIn.BrPP;                       % Br used for postprocessing [T]
    per.temphous = dataIn.HousingTemp;            % Housing Temperature [C]
    per.tempcuest = dataIn.EstimatedCopperTemp;   % Estimated Copper Temperatue [C]
    % torque and ripple penalization
    per.min_exp_torque = dataIn.MinExpTorque;      % minimum expected torque [Nm]
    per.max_exp_ripple = dataIn.MaxRippleTorque;    % maximum expected torque ripple in pu during opimization
    per.max_Cu_mass = dataIn.MaxCuMass;         % maximum expected copper mass [kg]
    if dataIn.LossEvaluationCheck == 1
        per.EvalSpeed = dataIn.EvalSpeed;
    end
    % a penalizing coefficient will be applied during the optimization to all
    % the machines with a lower torque or a higher torque ripple
    
    geo.BLKLABELSmaterials = {
        'Air';
        'Air';
        dataIn.SlotMaterial;
        dataIn.StatorMaterial;
        dataIn.RotorMaterial;
        dataIn.FluxBarrierMaterial;
        dataIn.ShaftMaterial;
        dataIn.RotorCondMaterial};
    
    geo.pont0 = dataIn.MinMechTol;  % thickness of the structural bridges at the airgap [mm]
    
    % Geometry
    geo.RotType = dataIn.TypeOfRotor;
    % 'Circular' is the Circular barrier type of rotor, for any number of barriers
    % 'ISeg' draws a rotor with the external I-shaped barrier and other
    %        Segmented (U-shaped) barriers
    % 'Seg' draws all segmented barriers
    % 'Fluid' draws barriers shaped according to fluid mechanics
    % 'SPM'     : Surface mounted permanent magnet motor
    geo.RemoveTMPfile = dataIn.RMVTmp;      % 'ON' of 'OFF' for remuving the motor folders in tmp
    
    geo.RaccBarrier='OFF';
    geo.DTrasl=0;
    
    geo.p  = dataIn.NumOfPolePairs;   % pole pairs
    geo.R  = dataIn.StatorOuterRadius;% stator outer radius [mm]
    geo.r  = dataIn.AirGapRadius;     % airgap radius [mm]
    geo.Ar = dataIn.ShaftRadius;      % shaft radius [mm]
    geo.g  = dataIn.AirGapThickness;  % airgap [mm]
    geo.l  = dataIn.StackLength;      % stack length [mm]
    
    % stator
    geo.q   = dataIn.NumOfSlots;      % stator slots per pole per phase
    if dataIn.SlotLayerPosCheck      % stator slot layer position, two solution are possible, 'over_under', and 'side_by_side'
        geo.slot_layer_pos = 'side_by_side';
    else
        geo.slot_layer_pos = 'over_under';
    end
    geo.lt  = dataIn.ToothLength;    % tooth length [mm]
    geo.acs = dataIn.StatorSlotOpen; % stator slot opening [p.u.]
    geo.wt  = dataIn.ToothWidth;       % tooth width [mm]
    
    if not(isfield(dataIn,'BarFillFac'))
        dataIn.BarFillFac = 1;
    end
    geo.BarFillFac=dataIn.BarFillFac;               % barrier filling factor
    
    geo.ttd = dataIn.ToothTangDepth;   % tooth tang depth [mm]
    geo.tta = dataIn.ToothTangAngle;   % tooth tang angle (mech degree)
    geo.SFR = dataIn.FilletCorner;     % fillet at the back corner of the slot [mm]
    
    % rotor
    if strcmp(geo.RotType,'SPM')
        geo.nlay = 1;
    else
        geo.nlay  = dataIn.NumOfLayers;    % number of layers
    end
    
    geo.racc_pont = 1 * geo.pont0;                  % radius of the fillet at the sides of inner bridges (if any) [mm]
    geo.ang_pont0 = geo.pont0 / geo.r * 180/pi;    % span of the arc corresponding to pont0 at the airgap radius [deg]
    geo.hfe_min = 2*geo.pont0;                      % min tickness of each steel flux guide
    % winding description
    geo.kcu = dataIn.SlotFillFactor;                % slot filling factor (net copper/slot area)
    tmp_avv=dataIn.WinMatr;
    geo.avv=[tmp_avv(:,1),[tmp_avv(2,2:end);tmp_avv(1,2:end)]];
    %     geo.avv = dataIn.WinMatr;
    
    geo.kracc = dataIn.PitchShortFac;       % pitch shortening factor (for end connections length estimation)
    geo.Ns = dataIn.TurnsInSeries;          % turns in series per phase (entire motor, one way, all poles in series)
    geo.ns = geo.q*6;                       % number of slot per pole pair
    geo.Nbob = geo.Ns/geo.p/(geo.ns/6)/size(geo.avv,1);  % conductors in slot per label
    if isfield(dataIn,'Qs')
        geo.Qs = dataIn.Qs;                     % number of stator slots in the FEMM simulation
    end
    geo.nmax = dataIn.OverSpeed; % overspeed [rpm]
    
    % permanent magnets
%     geo.Br = dataIn.Br;
%     geo.Hc = 1/(4e-7*pi)*geo.Br;
%     geo.rhoPM = 7600;   % Permanent Magnet mass density [Kg/m3] (for NdFeB magnets)
%     geo.Br_commercial=1.25; % [T] remanence of commercial NdFeB permanent magnets
    
    geo.lm = dataIn.ThicknessOfPM;
    geo.phi = dataIn.AngleSpanOfPM;
    % direction of magnetization in PMs of SPM with multiple segments of PM
    geo.PMdir = 'p';    % parallel direction
    geo.PMdir = 'r';    % radial direction
    
    % Mesh ratio (all_motor/air-gap)
    geo.K_mesh_MOOA = dataIn.Mesh_MOOA;    % optimization
    geo.K_mesh = dataIn.Mesh;         % post-processing and manual design
    
    % number of simulated rotor positions
    % MOOA means during the optimization
    geo.nsim_MOOA = dataIn.SimPoMOOA;         % simulated positions (6-1)
    geo.randFactor = dataIn.randFactor;       % Noise factor for position number reduction
    geo.delta_sim_MOOA = dataIn.RotPoMOOA;    % rotor position span [elt degrees]
    % evalx means the re-evaluation stage, with a finer resolution
    geo.nsim_singt = dataIn.SimPoFine;        % simulated positions (16-1)
    geo.delta_sim_singt = dataIn.RotPoFine;   % rotor position span [elt degrees]
    
    geo.x0=[0 0]; % geo.alpha = 0;
    geo.x0 = geo.r/cos(pi/2/geo.p);
    geo.dalpha_pu = dataIn.ALPHApu;
    geo.hc_pu = dataIn.HCpu;
    geo.dx = dataIn.DepthOfBarrier;
    
    geo.radial_ribs_eval = dataIn.RadRibCheck;  % if 1, the radial ribs are evaluated, else the radial ribs are setted in the GUI
    geo.pont = dataIn.RadRibEdit;
    %% bounds: limits of the search space
    % dalpha1 [p.u.]
    bounds_dalpha_1 = [dataIn.Alpha1Bou(1) dataIn.Alpha1Bou(2)];   % first angle [deg]
    % dalpha(2:nlay) [p.u.]
    if geo.nlay == 1
        bounds_dalpha = [dataIn.DeltaAlphaBou(1) dataIn.DeltaAlphaBou(2)]; % other angles [p.u.]
    else
        bounds_dalpha = ones(geo.nlay-1,1) * [dataIn.DeltaAlphaBou(1) dataIn.DeltaAlphaBou(2)]; % other angles [p.u.]
    end
    
    if strcmp(geo.RotType, 'SPM')
        % thickness of PM
        bounds_hc = geo.lm * [dataIn.hcBou(1) dataIn.hcBou(2)];
        geo.hc_pu = mean(bounds_hc);
    else
        % barrier ticknesses [p.u.]
        bounds_hc = ones(geo.nlay,1) * [dataIn.hcBou(1)  dataIn.hcBou(2)];
    end
    % barriers radial offset [p.u.]
    bounds_dx = ones(geo.nlay,1) * [dataIn.DfeBou(1) dataIn.DfeBou(2)];
    % remanence of the PMs
    bounds_Br = ones(geo.nlay,1) * [dataIn.BrBou(1) dataIn.BrBou(2)];
    % airgap
    bounds_g = [dataIn.GapBou(1)  dataIn.GapBou(2)];
    % Rotor radius
    bounds_xr = [dataIn.GapRadiusBou(1) dataIn.GapRadiusBou(2)];
    % Tooth width [mm]
    bounds_wt = [dataIn.ToothWiBou(1) dataIn.ToothWiBou(2)];
    % Tooth length [mm]
    bounds_lt = [dataIn.ToothLeBou(1) dataIn.ToothLeBou(2)];
    % stator slot opening [p.u.]
    bounds_acs = [dataIn.StatorSlotOpenBou(1) dataIn.StatorSlotOpenBou(2)];
    % tooth tang depth [mm]
    bounds_ttd = [dataIn.ToothTangDepthBou(1) dataIn.ToothTangDepthBou(2)];
    % phase angle of the current vector
    bounds_gamma = [dataIn.PhaseAngleCurrBou(1) dataIn.PhaseAngleCurrBou(2)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SETTINGS OF MODE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% RQ and bounds
    % RQ: vector of the n optimization inputs
    % bounds: n x 2 vector containing the boundaries of each input
    yes_vector = ones(geo.nlay,1);
    no_vector  = zeros(geo.nlay,1);
    if (strcmp(geo.RotType,'Fluid') || strcmp(geo.RotType,'Seg') || strcmp(geo.RotType,'Circular'))
        flag_dx = yes_vector*dataIn.DxBouCheck;
    else
        flag_dx = no_vector*dataIn.DxBouCheck;
    end
    
    
    if geo.nlay == 1
        bounds = [
            bounds_dalpha_1 dataIn.Dalpha1BouCheck
            bounds_dalpha   yes_vector(1)*dataIn.DalphaBouCheck
            bounds_hc       dataIn.hcBouCheck
            bounds_dx       flag_dx
            bounds_Br       dataIn.BrBouCheck
            bounds_g        dataIn.GapBouCheck
            bounds_xr       dataIn.AirgapRadiusBouCheck
            bounds_wt       dataIn.ToothWidthBouCheck
            bounds_lt       dataIn.ToothLengthBouCheck
            bounds_acs      dataIn.StatorSlotOpenBouCheck
            bounds_ttd      dataIn.ToothTangDepthBouCheck
            bounds_gamma    dataIn.GammaBouCheck];
    else
        bounds = [
            bounds_dalpha_1 dataIn.Dalpha1BouCheck
            bounds_dalpha   yes_vector(1:end-1)*dataIn.DalphaBouCheck
            bounds_hc       yes_vector*dataIn.hcBouCheck
            bounds_dx       flag_dx
            bounds_Br       yes_vector*dataIn.BrBouCheck
            bounds_g        dataIn.GapBouCheck
            bounds_xr       dataIn.AirgapRadiusBouCheck
            bounds_wt       dataIn.ToothWidthBouCheck
            bounds_lt       dataIn.ToothLengthBouCheck
            bounds_acs      dataIn.StatorSlotOpenBouCheck
            bounds_ttd      dataIn.ToothTangDepthBouCheck
            bounds_gamma    dataIn.GammaBouCheck];
    end
    
    filt_bounds = (bounds(:,3)==1);
    if geo.nlay == 1
        filt_bounds(2) = [];
        bounds(2,:) = [];
    end
    bounds = bounds(filt_bounds,1:2);
    
    %% OBJECTIVES
    objs = [per.min_exp_torque  dataIn.TorqueOptCheck
        per.max_exp_ripple      dataIn.TorRipOptCheck
        per.max_Cu_mass         dataIn.MassCuOptCheck];  % copper volume minimization
    filt_objs = (objs(:,2)==1);
    objs = objs(objs(:,2)==1,:);
    % eliminate unnecessary objs
    %objs = objs(objs(:,2)==1);
    per.objs=objs;
    
end

% names of the variables in RQ
RQnames{1} = 'dalpha';
if geo.nlay > 1
    for k = 2:geo.nlay
        RQnames{k} = 'dalpha';
    end
    kend = k;
    for k = kend+1:kend+geo.nlay
        RQnames{k} = 'hc';
    end
    kend = k;
    for k = kend+1:kend+geo.nlay
        RQnames{k} = 'dx';
    end
    kend = k;
    for k = kend+1:kend+geo.nlay
        RQnames{k} = 'Br';
    end
else
    RQnames{2} = 'hc';
    RQnames{3} = 'dx';
    RQnames{4} = 'Br';
    k = 4;
end
RQnames{k+1} = 'g';       % airgap
RQnames{k+2} = 'r';       % rotor radius
RQnames{k+3} = 'wt';      % tooth width
RQnames{k+4} = 'lt';      % tooth length
RQnames{k+5} = 'acs';     % stator slot opening [p.u.]
RQnames{k+6} = 'ttd';     % tooth tang depth [mm]
RQnames{k+7} = 'gamma';   % idq current phase angle

% eliminate unnecessary RQnames
RQnames = RQnames(filt_bounds);
geo.RQnames = RQnames;

% names of the MODE objectives 
OBJnames{1} = 'Torque';
OBJnames{2} = 'TorRip';
OBJnames{3} = 'MassCu';

% eliminate unnecessary OBJnames
OBJnames = OBJnames(filt_objs);
geo.OBJnames = OBJnames;

% Machine periodicity
Q = geo.ns*geo.p;                    % number of slots
geo.t = gcd(round(geo.ns*geo.p),geo.p);  % number of periods
% if ((6*geo.t/Q)>1)
%     % periodic machine
%     psCalc = 2*geo.p/geo.t;              % ps = # of poles of FEMM model
%     QsCalc = Q/geo.t;                    % Qs = # of slots of FEMM model
%     geo.tempWinTable = geo.avv;
%     geo.periodicity = 4;
% else
%     % anti-periodic machine
%     psCalc = geo.p/geo.t;                % ps = # of poles of FEMM model
%     QsCalc = Q/2/geo.t;                  % Qs = # of slots of FEMM model
%     geo.tempWinTable = [geo.avv -geo.avv];
%     geo.periodicity = 5;
% end

t2=gcd(round(geo.ns*geo.p),2*geo.p);
QsCalc=Q/t2;
psCalc=2*geo.p/t2;
if rem(psCalc,2)==0
    %periodic machine
    geo.tempWinTable = geo.avv;
    geo.periodicity = 4;
else
    %anti-periodic machine
    geo.tempWinTable = [geo.avv -geo.avv];
    geo.periodicity = 5;
end
if isfield(geo,'Qs')  % Qs set in the GUI
    geo.ps = round(psCalc*geo.Qs/QsCalc);
    % Check for periodicity: if ps is odd, the simulated portion of motor
    % is anti-periodic; if ps is even, the simulated portion is periodic.
    if rem(geo.ps,2)==0
        geo.periodicity = 4;   % periodic configuration
    else
        geo.periodicity = 5;   % anti-periodic configuration
    end
    geo.tempWinTable = geo.avv;
else
    geo.Qs = QsCalc;
    geo.ps = psCalc;
end
    



% Yield strength definition for rotor laminations
% if strcmp(char(geo.BLKLABELSmaterials(5)), 'vacodur-opt-mec')>0
%     geo.sigma_max = 390;    % [MPa]
%     geo.rhoFE = 8120;       % mass density [Kg/m3]
% elseif strcmp(char(geo.BLKLABELSmaterials(5)), 'Transil270-35')>0
%     geo.sigma_max = 180;    % [MPa]
%     geo.rhoFE = 7800;       % mass density [Kg/m3]
% elseif strcmp(char(geo.BLKLABELSmaterials(5)), '10JNEX900')>0
%     geo.sigma_max = 604;    % [MPa]
%     geo.rhoFE = 7490;       % mass density [Kg/m3]
%     geo.loss.alpha = 1.15746;    % coefficient of Steinmetz equation
%     geo.loss.beta = 1.78672;     % coefficient of Steinmetz equation
%     geo.loss.kh = 0.00571782;    % coefficient of Steinmetz equation
%     geo.loss.ke = 3.71539*1e-006;% coefficient of Steinmetz equation
% elseif strcmp(char(geo.BLKLABELSmaterials(5)), 'M250-35A')>0
%     geo.sigma_max = 455;    % [MPa]
%     geo.rhoFE = 7600;       % mass density [Kg/m3]
%     geo.loss.alpha = 1.23089;    % coefficient of Steinmetz equation
%     geo.loss.beta = 1.79026;     % coefficient of Steinmetz equation
%     geo.loss.kh = 0.00777985;    % coefficient of Steinmetz equation
%     geo.loss.ke = 3.14545e-005;% coefficient of Steinmetz equation
% elseif strcmp(char(geo.BLKLABELSmaterials(5)), 'M530-65A-OK')>0
%     geo.sigma_max = 285;    % [MPa]
%     geo.rhoFE = 7700;       % mass density [Kg/m3]
%     geo.loss.alpha = 1;    % coefficient of Steinmetz equation
%     geo.loss.beta = 1.7265;     % coefficient of Steinmetz equation
%     geo.loss.kh = 0.0413163;    % coefficient of Steinmetz equation
%     geo.loss.ke = 0;% coefficient of Steinmetz equation
% elseif strcmp(char(geo.BLKLABELSmaterials(5)), 'Hyperco 0.006') % added for support
%     geo.sigma_max = 717;    % [MPa]
%     geo.rhoFE = 8120;       % mass density [kg/m3]
%     geo.loss.alpha = 1.15293;    % coefficient of Steinmetz equation
%     geo.loss.beta = 1.7224;      % coefficient of Steinmetz equation
%     geo.loss.kh = 0.00737903;    % coefficient of Steinmetz equation
%     geo.loss.ke = 9.26301e-006;  % coefficient of Steinmetz equation
% else
%     geo.sigma_max = 200;    % [MPa]
%     geo.rhoFE = 7800;       % mass density [Kg/m3]
%     geo.loss.alpha = 0;    % coefficient of Steinmetz equation
%     geo.loss.beta = 0;     % coefficient of Steinmetz equation
%     geo.loss.kh = 0;    % coefficient of Steinmetz equation
%     geo.loss.ke = 0;% coefficient of Steinmetz equation
%     %     disp('the rotor material Yield strength was undefined; it has been set to 200 MPa')
% end

mat = assign_mat_prop(dataIn);

