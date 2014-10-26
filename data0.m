function [bounds, geo, per] = data0(dataIn)
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
%% data0.m:
% manual input data used as default by the graphic input GUI

%% the data are organized by data structures
% per: performance data
% geo: geometrical data

if (nargin<1)
    
    %% MANUAL SETTINGS
    %% per: performance target
    per.Loss = 65;             % admitted Joule loss [W]
    per.tempcu = 100;          % Copper Temperature [C]
    per.overload = 2;          % current overload factor used for optimization (1 means Joule loss = per.Loss)
    
    % torque and ripple penalization
    % a penalizing coefficient will be applied during the optimization to all
    % the machines with a lower torque or a higher torque ripple
    per.min_exp_torque=  10;   % minimum expected torque [Nm]
    per.max_exp_ripple = 8;    % maximum expected torque ripple in pu during opimization
    
    %% geo: geometry
    geo.BLKLABELSmaterials={
        'Air';
        'Air';
        'Copper';
        'M530_65A_MagNet';
        'M530_65A_MagNet';
        'Bonded-Magnet';
        'M530_65A_MagNet';
        'Copper'};
    
    geo.pont0 = 0.4;    % thickness of the structural bridges at the airgap
    % AND minimum mechanical tolerance [mm]
    
    % Type of barriers
    % 'Circular': Circular barriers
    % 'ISeg'    : I-shaped exterior barrier plus Seg(mented) or U-shaped barriers
    % 'Seg'     : All Seg(mented) barriers
    % 'Fluid'   : barriers shaped according to fluid flow lines
    geo.RotType = 'Circular';
    
    geo.RemoveTMPfile = 'ON';      % 'ON' of 'OFF' for remuving the motor folders in tmp
    
    geo.p  = 2;         % pole pairs
    geo.R  = 67.5;      % stator outer radius [mm]
    geo.r  = 0.6*geo.R-0.275;  % airgap radius [mm]
    geo.Ar = geo.r*(sqrt(2)-1);     % shaft radius [mm]
    geo.g  = 0.275;     % airgap [mm]
    geo.l  = 101;       % stack length [mm]
    
    % stator
    geo.q   = 2;        % stator slots per pole per phase
    
    % tooth length and width
    b  = 0.5;       % Bgap/Bfe,yoke  (back-iron p.u. size)
    bt = 0.5;       % Bgap/Bfe,tooth (tooth p.u. size)
    geo.lt = geo.R * (1-geo.r/geo.R*(1+b/geo.p));
    geo.wt = pi/(3*geo.q*geo.p)*geo.r*bt;
    
    geo.acs  = 0.25;    % stator slot opening [p.u.]
    geo.ttd = 0.75;     % tooth tang depth [mm]
    geo.tta = 25;       % tooth tang angle (mech degree)
    geo.SFR = 1;        % fillet at the back corner of the slot [mm]
    
    % rotor
    geo.nlay  = 3;      % number of layers
    
    geo.racc_pont = 1 * geo.pont0;                  % radius of the fillet at the sides of inner bridges (if any) [mm]
    geo.ang_pont0 = geo.pont0 / geo.r * 180/pi;     % span of the arc corresponding to pont0 at the airgap radius [deg]
    geo.hfe_min = 2*geo.pont0;                      % minimum tickness of each steel flux guide
    % winding description
    geo.kcu = 0.407;                                % slot filling factor (net copper/slot area)
    geo.avv   = [1 1 -3 -3 2 2                      % winding scheme
                 1 1 -3 -3 2 2];
    
    geo.kracc = 9/9;                                % pitch shortening factor (for end connections length estimation)
    geo.Ns    = 122;                                % turns in series per phase (entire motor, one way, all poles in series)
    geo.ns = geo.q*6;                               % number of slot per pole pair
    geo.Nbob  = geo.Ns/geo.p/(geo.ns/6)/size(geo.avv,1);  % conductors in slot per label
    geo.nmax = 8000;                                % overspeed [rpm]
    
    % permanent magnets (default: Br = 0)
    geo.Br = 0.0;
    geo.Hc = 1/(4e-7*pi)*geo.Br;
    geo.rhoPM = 7600;                               % Permanent Magnet mass density [Kg/m3] (for NdFeB magnets)
    geo.Br_commercial=1.25;                         % [T] remanence of commercial NdFeB permanent magnets
    
    % Settings that vary if MOOA or not (MOOA means during the
    % optimization)
    % Mesh ratio (all_motor/air-gap)
    geo.K_mesh_MOOA = 4;    % optimization
    geo.K_mesh = 4;         % post-processing and manual design
    
    % number of simulated rotor positions
    geo.nsim_MOOA = 5;          % simulated positions
    geo.delta_sim_MOOA = 360/(geo.q*6);    % rotor position span [elt degrees]
    % evalx means the re-evaluation stage, with a finer resolution
    geo.nsim_singt = 30;        % simulated positions
    geo.delta_sim_singt = 60;   % rotor position span [elt degrees]
    
    % rotor geometry (subject to optimization: default set to zero)
    geo.x0=[0 0]; geo.alpha = 0;
    geo.dalpha_pu = [0.45 0.22*ones(1,geo.nlay-1)];
    geo.hc_pu = 0.5*ones(1,geo.nlay);
    geo.dx = ones(1,geo.nlay) * 0;
    
    %% bounds: limits of the search space
    % dalpha1 [p.u.]
    bounds_dalpha_1 = [22.5 45]/90;   % first angle [deg]
    % dalpha(2:nlay) [p.u.]
    bounds_dalpha = ones(geo.nlay-1,1) * [0.5/geo.nlay 0.5];% other angles [p.u.]
    % barrier ticknesses [p.u.]
    bounds_hc = ones(geo.nlay,1) * [0.2 1];
    % barriers radial offset [p.u.]
    bounds_dx=ones(geo.nlay,1) * [-0.75 0.75];
    % remanence of the PMs
    bounds_Br = [0.2 0.8];
    % airgap
    bounds_g = geo.g * [0.7  1.5];
    % Rotor radius
    bounds_xr = geo.r * [0.8  1.2];
    % Tooth width [mm]
    bounds_wt = geo.wt * [0.75  1.3];
    % Tooth length [mm]
    bounds_lt = geo.lt * [0.8 1.2];
    % phase angle of the current vector
    bounds_gamma = [40 75];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% RQ and bounds: define here the variables to be optimized
    %% RQ is the vector of the optimization inputs (variable size n)
    %% bounds is the n x 2 vector containig the boundaries of each input
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yes_vector = ones(geo.nlay,1);
    no_vector = zeros(geo.nlay,1);
    if (strcmp(geo.RotType,'Fluid') || strcmp(geo.RotType,'Seg'))
        flag_dx = yes_vector;
    else
        flag_dx = no_vector;
    end
    
    % Names of the variables in RQ
    names{1} = 'dalpha_1';
    if geo.nlay > 1
        for k = 2:geo.nlay
            names{k} = 'dalpha';
        end
        kend = k;
    else
        kend = 2;
    end
    for k = kend+1:kend+geo.nlay
        names{k} = 'hc';
    end
    kend = k;
    for k = kend+1:kend+geo.nlay
        names{k} = 'dx';
    end
    names{k+1} = 'Br';      % magnets
    names{k+2} = 'g';      % airgap
    names{k+3} = 'r';      % rotor radius
    names{k+4} = 'wt';      % tooth width
    names{k+5} = 'lt';      % tooth length
    names{k+6} = 'gamma';   % idq current phase angle
    
    % bounds
    % put 1 or 0 to add/remove variables from RQ and bounds
    
    % example:
    bounds = [
        bounds_dalpha_1 1
        bounds_dalpha   yes_vector(1:end-1)
        bounds_hc       yes_vector
        bounds_dx       flag_dx
        bounds_Br       0
        bounds_g        0
        bounds_xr       0
        bounds_wt       0
        bounds_lt       0
        bounds_gamma    1];
    
    % eliminate unnecessary rows
    filt_bounds = (bounds(:,3)==1);
    bounds = bounds(filt_bounds,1:2);
    % eliminate the unnecessary names
    names = names(filt_bounds);
    
    geo.RQnames = names;
    
else
    
    %% READ INPUTS FROM THE GUI
    
    % main performance target
    per.Loss = dataIn.AdmiJouleLosses;             % admitted Joule loss [W]
    per.tempcu = dataIn.CopperTemp;           % Copper Temperature [C]
    %     per.Vdc = dataIn.DCVoltage;              % dc-link voltage [V]
    per.overload = dataIn.CurrOverLoad;           % current overload factor used for optimization (1 means Joule loss = per.Loss)
    
    % torque and ripple penalization
    per.min_exp_torque = dataIn.MinExpTorque;      % minimum expected torque [Nm]
    per.max_exp_ripple = dataIn.MaxRippleTorque;    % maximum expected torque ripple in pu during opimization
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
    geo.RotType = dataIn.TypeOfRotor;       % 'Circular' is the Circular barrier type of rotor, for any number of barriers
    % 'ISeg' draws a rotor with the external I-shaped barrier and other
    %        Segmented (U-shaped) barriers
    % 'Seg' draws all segmented barriers
    % 'Fluid' draws barriers shaped according to fluid mechanics
    geo.RemoveTMPfile = dataIn.RMVTmp;      % 'ON' of 'OFF' for remuving the motor folders in tmp
    
    geo.p  = dataIn.NumOfPolePairs;   % pole pairs
    geo.R  = dataIn.StatorOuterRadius;% stator outer radius [mm]
    geo.r = dataIn.AirGapRadius;     % airgap radius [mm]
    geo.Ar = dataIn.ShaftRadius;      % shaft radius [mm]
    geo.g  = dataIn.AirGapThickness;  % airgap [mm]
    geo.l  = dataIn.StackLength;      % stack length [mm]
    
    % stator
    geo.q   = dataIn.NumOfSlots;      % stator slots per pole per phase
    geo.lt  = dataIn.ToothLength;    % tooth length [mm]
    geo.acs = dataIn.StatorSlotOpen; % stator slot opening [p.u.]
    geo.wt  = dataIn.ToothWidth;       % tooth width [mm]
    
    geo.ttd = dataIn.ToothTangDepth;   % tooth tang depth [mm]
    geo.tta = dataIn.ToothTangAngle;   % tooth tang angle (mech degree)
    geo.SFR = dataIn.FilletCorner;     % fillet at the back corner of the slot [mm]
    
    % rotor
    geo.nlay  = dataIn.NumOfLayers;    % number of layers
    
    geo.racc_pont = 1 * geo.pont0;                  % radius of the fillet at the sides of inner bridges (if any) [mm]
    geo.ang_pont0 = geo.pont0 / geo.r * 180/pi;    % span of the arc corresponding to pont0 at the airgap radius [deg]
    geo.hfe_min = 2*geo.pont0;                      % min tickness of each steel flux guide
    % winding description
    geo.kcu = dataIn.SlotFillFactor;                % slot filling factor (net copper/slot area)
    geo.avv = dataIn.WinMatr;
    
    geo.kracc = dataIn.PitchShortFac;       % pitch shortening factor (for end connections length estimation)
    geo.Ns = dataIn.TurnsInSeries;          % turns in series per phase (entire motor, one way, all poles in series)
    geo.ns = geo.q*6;                       % number of slot per pole pair
    geo.Nbob = geo.Ns/geo.p/(geo.ns/6)/size(geo.avv,1);  % conductors in slot per label
    geo.nmax = dataIn.OverSpeed; % overspeed [rpm]
    
    % permanent magnets
    geo.Br = dataIn.Br;
    geo.Hc = 1/(4e-7*pi)*geo.Br;
    geo.rhoPM = 7600;   % Permanent Magnet mass density [Kg/m3] (for NdFeB magnets)
    geo.Br_commercial=1.25; % [T] remanence of commercial NdFeB permanent magnets
    
    % Mesh ratio (all_motor/air-gap)
    geo.K_mesh_MOOA = dataIn.Mesh_MOOA;    % optimization
    geo.K_mesh = dataIn.Mesh;         % post-processing and manual design
    
    % number of simulated rotor positions
    % MOOA means during the optimization
    geo.nsim_MOOA = dataIn.SimPoMOOA;          % simulated positions (6-1)
    geo.delta_sim_MOOA = dataIn.RotPoMOOA;    % rotor position span [elt degrees]
    % evalx means the re-evaluation stage, with a finer resolution
    geo.nsim_singt = dataIn.SimPoFine;        % simulated positions (16-1)
    geo.delta_sim_singt = dataIn.RotPoFine;   % rotor position span [elt degrees]
    
    geo.x0=[0 0]; geo.alpha = 0;
    geo.dalpha_pu = [0.45 0.22*ones(1,geo.nlay-1)];
    geo.hc_pu = 0.5*ones(1,geo.nlay);
    geo.dx = ones(1,geo.nlay) * 0;
    
    %% bounds: limits of the search space
    % dalpha1 [p.u.]
    bounds_dalpha_1 = [dataIn.Alpha1Bou(1) dataIn.Alpha1Bou(2)];   % first angle [deg]
    % dalpha(2:nlay) [p.u.]
    bounds_dalpha = ones(geo.nlay-1,1) * [dataIn.DeltaAlphaBou(1) dataIn.DeltaAlphaBou(2)]; % other angles [p.u.]
    % barrier ticknesses [p.u.]
    bounds_hc = ones(geo.nlay,1) * [dataIn.hcBou(1)  dataIn.hcBou(2)];
    % barriers radial offset [p.u.]
    bounds_dx = ones(geo.nlay,1) * [dataIn.DfeBou(1) dataIn.DfeBou(2)];
    % remanence of the PMs
    bounds_Br = [dataIn.BrBou(1) dataIn.BrBou(2)];
    % airgap
    bounds_g = [dataIn.GapBou(1)  dataIn.GapBou(2)];
    % Rotor radius
    bounds_xr = [dataIn.GapRadiusBou(1) dataIn.GapRadiusBou(2)];
    % Tooth width [mm]
    bounds_wt = [dataIn.ToothWiBou(1) dataIn.ToothWiBou(2)];
    % Tooth length [mm]
    bounds_lt = [dataIn.ToothLeBou(1) dataIn.ToothLeBou(2)];
    % phase angle of the current vector
    bounds_gamma = [dataIn.PhaseAngleCurrBou(1) dataIn.PhaseAngleCurrBou(2)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% RQ and bounds: define here the variables to be optimized
    %% RQ is the vector of the optimization inputs (variable size n)
    %% bounds is the n x 2 vector containig the boundaries of each input
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yes_vector = ones(geo.nlay,1);
    no_vector = zeros(geo.nlay,1);
    if (strcmp(geo.RotType,'Fluid') || strcmp(geo.RotType,'Seg'))
        flag_dx = yes_vector*dataIn.DxBouCheck;
    else
        flag_dx = no_vector*dataIn.DxBouCheck;
    end
    
    % names of the variables in RQ
    names{1} = 'dalpha_1';
    if geo.nlay > 1
        for k = 2:geo.nlay
            names{k} = 'dalpha';
        end
        kend = k;
    else
        kend = 2;
    end
    for k = kend+1:kend+geo.nlay
        names{k} = 'hc';
    end
    kend = k;
    for k = kend+1:kend+geo.nlay
        names{k} = 'dx';
    end
    names{k+1} = 'Br';      % magnets
    names{k+2} = 'g';      % airgap
    names{k+3} = 'r';      % rotor radius
    names{k+4} = 'wt';      % tooth width
    names{k+5} = 'lt';      % tooth length
    names{k+6} = 'gamma';   % idq current phase angle
    
    bounds = [
        bounds_dalpha_1 dataIn.Dalpha1BouCheck
        bounds_dalpha   yes_vector(1:end-1)*dataIn.DalphaBouCheck
        bounds_hc       yes_vector*dataIn.hcBouCheck
        bounds_dx       flag_dx
        bounds_Br       dataIn.BrBouCheck
        bounds_g        dataIn.GapBouCheck
        bounds_xr       dataIn.AirgapRadiusBouCheck
        bounds_wt       dataIn.ToothWidthBouCheck
        bounds_lt       dataIn.ToothLengthBouCheck
        bounds_gamma    dataIn.GammaBouCheck];
    
    filt_bounds = (bounds(:,3)==1);
    bounds = bounds(filt_bounds,1:2);
    % eliminate the unnecessary names
    names = names(filt_bounds);
    geo.RQnames = names;
    
end

% Yield strength definition for rotor laminations
if strcmp(char(geo.BLKLABELSmaterials(5)), 'vacodur-opt-mec')>0
    geo.sigma_max = 390;    % [MPa]
    geo.rhoFE = 8120;       % mass density [Kg/m3]
elseif strcmp(char(geo.BLKLABELSmaterials(5)), 'Transil270-35')>0
    geo.sigma_max = 180;    % [MPa]
    geo.rhoFE = 7800;       % mass density [Kg/m3]
elseif strcmp(char(geo.BLKLABELSmaterials(5)), '10JNEX900')>0
    geo.sigma_max = 604;    % [MPa]
    geo.rhoFE = 7490;       % mass density [Kg/m3]
else
    geo.sigma_max = 200;    % [MPa]
    geo.rhoFE = 7800;       % mass density [Kg/m3]
    %     disp('the rotor material Yield strength was undefined; it has been set to 200 MPa')
end

