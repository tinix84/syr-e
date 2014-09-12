function [bounds, geo, per] = data1(dataIn)
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
%% data0.m: file containing the input data

%% the data are organized by data structures 
% per: performance data
% mat: materials names
% BLKLABELS: refers to the blocks of the FEMM model
% geo: geometrical data

% main performance target
per.Loss = dataIn.AdmiJouleLosses;             % admitted Joule loss [W]
per.tempcu = dataIn.CopperTemp;           % Copper Temperature [C]
per.Vdc = dataIn.DCVoltage;              % dc-link voltage [V]
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

% Yield strength definition for rotor laminations
if strcmp(char(geo.BLKLABELSmaterials(5)), 'vacodur-opt-mec')>0
    geo.sigma_max = 390;    % [MPa]
elseif strcmp(char(geo.BLKLABELSmaterials(5)), 'Transil270-35')>0
    geo.sigma_max = 180;    % [MPa]
elseif strcmp(char(geo.BLKLABELSmaterials(5)), '10JNEX900')>0
    geo.sigma_max = 604;    % [MPa]
else
    geo.sigma_max = 200;    % [MPa]
%     disp('the rotor material Yield strength was undefined; it has been set to 200 MPa')
end

geo.pont0 = dataIn.MinMechTol;  % thickness of the structural bridges at the airgap [mm]

% Mesh ratio (all_motor/air-gap)
geo.K_mesh = dataIn.Mesh;
geo.K_mesh_MOOA = dataIn.Mesh_MOOA;

% Geometry 
geo.RotType = dataIn.TypeOfRotor;       % 'Circular' is the Circular barrier type of rotor, for any number of barriers
 % 'ISeg' draws a rotor with the external I-shaped barrier and other 
    %        Segmented (U-shaped) barriers
    % 'Seg' draws all segmented barriers
    % 'Fluid' draws barriers shaped according to fluid mechanics
    
geo.RemoveTMPfile = dataIn.RMVTmp;      % 'ON' of 'OFF' for remuving the motor folders in tmp

geo.p  = dataIn.NumOfPolePairs;   % pole pairs
geo.r  = dataIn.StatorOuterRadius;% stator outer radius [mm]
geo.xr = dataIn.AirGapRadius;     % airgap radius [mm]
geo.Ar = dataIn.ShaftRadius;      % shaft radius [mm]
geo.g  = dataIn.AirGapThickness;  % airgap [mm]
geo.l  = dataIn.StackLength;      % stack length [mm]

% stator
geo.q   = dataIn.NumOfSlots;      % stator slots per pole per phase
geo.lt   = dataIn.ToothLength;    % tooth length [mm]
geo.acs  = dataIn.StatorSlotOpen; % stator slot opening [p.u.]
geo.b    = 0.659;           % Bgap/Byoke (determines tooth width^-1, yoke width^-1)
geo.kt   = 0.8; % tooth scaling factor kt = Byoke/Btooth (further determines tooth width^-1)
geo.wt = dataIn.ToothWidth; % tooth width [mm]

geo.ttd = dataIn.ToothTangDepth;   % tooth tang depth [mm]
geo.tta = dataIn.ToothTangAngle;   % tooth tang angle (mech degree)
geo.SFR = dataIn.FilletCorner;     % fillet at the back corner of the slot [mm]

% rotor
geo.nlay  = dataIn.NumOfLayers;    % number of layers

geo.racc_pont = 1 * geo.pont0;                  % radius of the fillet at the sides of inner bridges (if any) [mm] 
geo.ang_pont0 = geo.pont0 / geo.xr * 180/pi;    % span of the arc corresponding to pont0 at the airgap radius [deg]
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

%% number of simulated rotor positions
% MOOA means during the optimization
geo.nsim_MOOA = dataIn.SimPoMOOA;          % simulated positions (6-1)
geo.delta_sim_MOOA = dataIn.RotPoMOOA;    % rotor position span [elt degrees] 
% evalx means the re-evaluation stage, with a finer resolution
geo.nsim_singt = dataIn.SimPoFine;        % simulated positions (16-1)
geo.delta_sim_singt = dataIn.RotPoFine;   % rotor position span [elt degrees] 

geo.x0 = [0 0]; geo.alpha = 0;
geo.dalpha = zeros(1,geo.nlay+1);
geo.hc_pu = zeros(1,geo.nlay);
var.hc_pu = geo.hc_pu;

%% Limits of the search space
% dalpha and dalpha1 [deg] 
bounds_dalpha_1 = [dataIn.Alpha1Bou(1) dataIn.Alpha1Bou(2)];   % first angle [deg]
temp_dalpha =  1/(geo.nlay) * ones(1,geo.nlay);
bounds_dalpha = temp_dalpha' * [dataIn.DeltaAlphaBou(1) dataIn.DeltaAlphaBou(2)];                               % other angles [p.u.]                               
% barrier ticknesses [p.u.]
temp_hc =  ones(1,geo.nlay);    
bounds_hc = temp_hc' * [dataIn.hcBou(1)  dataIn.hcBou(2)];
%Spessore dei ferri (percentuale da togliere o aggiungere all'aria)
temp_Dfe=  ones(1,geo.nlay);
bounds_Dfe = temp_Dfe' * [dataIn.DfeBou(1) dataIn.DfeBou(2)];
% remanence of the PMs (not used for SyR optimization)
bounds_Br = [0.2 0.8];
% airgap radius
% x0 = geo.xr/geo.r;
% bounds_x = x0 * [0.8  1.15];
% Split ratio
% bounds_sr = [0.4 0.6];
% Tooth length [mm] 
% bounds_lt = [10 14];
% Tooth scaling factor kt
% bounds_kt = [0.55 0.95];
% phase angle of the current vector
bounds_gamma = [dataIn.PhaseAngleCurrBou(1) dataIn.PhaseAngleCurrBou(2)];
% Numero di barriere di flusso ('nlayers')
% bounds_nlay = [0.51 4.49];
% p.u. displacement of the 1st layer (Seg only)
% if (strcmp(geo.RotType,'Seg'))
%     bounds_Delta_X = [dataIn.DeltaXBou(1) dataIn.DeltaXBou(2)];
% end

%% bounds of the search space
if (strcmp(geo.RotType,'Fluid') || strcmp(geo.RotType,'3U_Dfe')||strcmp(geo.RotType,'Seg'))
    bounds = [
        bounds_dalpha_1
        bounds_dalpha
        bounds_hc
        bounds_Dfe
        bounds_gamma];
else
    bounds = [
        bounds_dalpha_1
        bounds_dalpha
        bounds_hc
        bounds_gamma];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%