%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% data0.m: file containing the input data

%% the data are organized by data structures 
% per: performance data
% mat: materials names
% BLKLABELS: refers to the blocks of the FEMM model
% geo: geometrical data

% main performance target
per.Loss = 192;             % admitted Joule loss [W]
per.tempcu = 100;           % Copper Temperature [C]
per.Vdc = 140;              % dc-link voltage [V]
per.overload = 2;           % current overload factor used for optimization (1 means Joule loss = per.Loss) 

% torque and ripple penalization
per.min_exp_torque= 4;      % minimum expected torque [Nm]
per.max_exp_ripple = 10;    % maximum expected torque ripple in pu during opimization
% a penalizing coefficient will be applied during the optimization to all
% the machines with a lower torque or a higher torque ripple

% block material organization
%BLKLABELS.materials = {air
%                       air
%                       slot material
%                       Stator yoke and tooth material
%                       Rotor yoke material
%                       flux barrier material
%                       shaft material
%                       rotor conductor winding material (not modeled)
%                       }



BLKLABELS.materials={
    'Air';
    'Air';
    'Copper';
    'S18-Corea';
    'Transil270-35';
    'Bonded-Magnet';
    'Transil270-35';
    'Copper'};

% Mesh ratio (all_motor/air-gap)
geo.K_mesh = 4;
geo.K_mesh_MOOA = 4;
geo.scalingfactor = 1;                % mm to FEMM units

% Geometry
geo.RotType='I2U';       % 3C is the Circulare barrier type of rotor, for any number of barriers

geo.p  = 5;                           % pole pairs
geo.r  = 50.5*geo.scalingfactor;      % stator outer radius [mm]
geo.xr = 29.29*geo.scalingfactor;     % airgap radius [mm]
geo.Ar = 10*geo.scalingfactor;         % shaft radius [mm]
geo.g  = 0.5;                         % airgap [mm]
geo.l  = 65;                          % stack length [mm]

% stator
geo.ns   = 2.4;                        % stator slots per pole pair
geo.lt   = 14.55*geo.scalingfactor;   % tooth length [mm]
geo.acs  = 0.3;                       % stator slot opening [p.u.]
geo.b    = 0.42;                      % Bgap/Btooth (determines tooth width^-1, yoke width^-1)
geo.kt   = 0.8;                       % tooth scaling factor kt = Byoke/Btooth (further determines tooth width^-1)

geo.ttd   = 0.76*geo.scalingfactor;    % tooth tang depth [mm]
geo.tta = 10;                       % tooth tang angle (mech degree)
geo.SFR = 1.5;                        % fillet at the back corner of the slot [mm]

% rotor
geo.nlay  = 3;                                 % number of layers
geo.pont0 = 0.4;                                % thickness of the structural bridges at the airgap [mm]
geo.racc_pont = 1 * geo.pont0;                  % radius of the fillet at the sides of inner bridges (if any) [mm] 
geo.ang_pont0 = geo.pont0 / geo.xr * 180/pi;    % span of the arc corresponding to pont0 at the airgap radius [deg]
geo.hfe_min=2*geo.pont0;                            % min tickness of each steel flux guide
% winding description
geo.kcu = 0.407;                      % slot filling factor (net copper/slot area)
geo.avv   = [ 1 1 -3 -3 2 2           % winding scheme
              1 1 -3 -3 2 2];
% geo.avv=[1,2,-2,-3,3,1,-1,-2,2,3,-3,-1;1,-1,-2,2,3,-3,-1,1,2,-2,-3,3];

geo.kracc = 6/6;                      % pitch shortening factor (for end connections length estimation)
geo.Ns    = 122;                       % turns in series per phase (entire motor, one way, all poles in series)

geo.Nbob  = geo.Ns/geo.p/(geo.ns/6)/size(geo.avv,1);  % conductors in slot per label
geo.nmax = 8000;                     % overspeed [rpm]

% permanent magnets
geo.Br = 0.0;
geo.Hc = 1/(4e-7*pi)*geo.Br;

%% number of simulated rotor positions
% MOOA means during the optimization
geo.nsim_MOOA = 6;          % simulated positions (6-1)
geo.delta_sim_MOOA = 30;    % rotor position span [elt degrees] 
% evalx means the re-evaluation stage, with a finer resolution
geo.nsim_singt = 31;        % simulated positions (16-1)
geo.delta_sim_singt = 60;   % rotor position span [elt degrees] 

geo.x0=[0 0]; geo.alpha = 0;
geo.dalpha = zeros(1,geo.nlay+1);
% geo.magpu = ones(1,geo.nlay);
geo.hc_pu = zeros(1,geo.nlay);
var.hc_pu = geo.hc_pu;
% geo = calc_alpha_hc_delta_x0_2(geo);

%% Limits of the search space

% dalpha and dalpha1 [deg] 
bounds_dalpha_1 = [180/(2*geo.p*geo.nlay) (360/(4*geo.p))/2*1.2];   % first angle [deg]
temp_dalpha =  1/(geo.nlay) * ones(1,geo.nlay);
bounds_dalpha = temp_dalpha' * [1 2];                               % other angles [p.u.]                               
% barrier ticknesses [p.u.]
temp_hc =  ones(1,geo.nlay);    
bounds_hc = temp_hc' * [0.2 1];
%Spessore dei ferri (percentuale da togliere o aggiungere all'aria)
temp_Dfe=  ones(1,geo.nlay);
bounds_Dfe=temp_Dfe' * [-0.75 0.75];
% remanence of the PMs (not used for SyR optimization)
bounds_Br = [0.2 0.8];
% airgap radius
x0 = geo.xr/geo.r;
bounds_x = x0 * [0.8  1.15];
% phase angle of the current vector
bounds_gamma = [20 80];
% Numero di barriere di flusso ('nlayers')
% bounds_nlay = [0.51 4.49];
% p.u. displacement of the 1st layer (3U only)
if (strcmp(geo.RotType,'3U'))
    bounds_Delta_X=[0 1];
end

%% bounds of the search space
if (strcmp(geo.RotType,'3U'))
    bounds = [
        bounds_dalpha_1
        bounds_dalpha
        bounds_hc
        bounds_Delta_X  % 3U only
        bounds_gamma];
elseif (strcmp(geo.RotType,'Fluid'))
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