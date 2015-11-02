%% VERSIONE 20 11 2011
% RunSimulation
% Input parameters in struct Sim:
%                               .n = 15000;          %rpm
%                               .Id_ = [ ... ];      % Id vector (1xn)
%                               .Iq_ = [ ... ];      % Iq vector (1xn)
%                               .sim_angle = 360;    % simulated angle (elt deg)
% Output:
% n x   \CaseBackup\MachineName_Id_Iq.mn      one MagNet Case per each (Id,Iq)pair
% n x   \CaseBackup\MachineName_Id_Iq.mat     Processed waveforms, single case
% 1 x   \Mat\map_max(Id)_max(Iq)_caso.mat     Fd,Fq,T, one point per simulation


% addpath M
switch_computer
load('temp_files\BATCH_SETUP.mat')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Costruzione dell'insieme di simulazioni da fare
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% La variabile (Sim) contiene TUTTE le indicazioni per TUTTE le
% simulazioni
% La variabile (Cas) contiene le specifiche della singola simulazione e'
% definita in "RunS_RunCurve.m"
Sim.temperature = In.temperature;
% Riconosce size(Imod) >< size(gamma)
if (length(In.Imod) > length(In.gamma))
    % (A)
    Sim.Imod  = In.Imod;
    Sim.gamma = In.gamma(1) * ones(size(In.Imod));
elseif length(In.Imod) < length(In.gamma)
    % (B)
    Sim.Imod  = In.Imod(1) * ones(size(In.gamma));
    Sim.gamma = In.gamma;
else
    % (C)
    Sim.Imod  = In.Imod;
    Sim.gamma = In.gamma;
end

% Calcola anche Id e Iq per ogni caso
Sim.Id = Sim.Imod .* cos(Sim.gamma*pi/180);
Sim.Iq = Sim.Imod .* sin(Sim.gamma*pi/180);

%% Impostazioni per lo skewing
Sim.skew       = In.skew;
Sim.skew_Nstep = In.skew_Nstep;
% Sim.skew_angle = In.skew_angle;  % Skew_New

%% Costruzione dei diversi casi in funzione dei parametri di (Cas)
% la variabile (In) non serve piu'.
if (Sim.skew)
    % Vettore skew_Imod
    Sim.skew_Imod = ones(Sim.skew_Nstep,1) * Sim.Imod;
    temp = In.skew_angle/Sim.skew_Nstep*(0:1:Sim.skew_Nstep-1);
    Sim.skew_angle = (temp - mean(temp))';   % Skew_New
    Sim.skew_gamma = ones(Sim.skew_Nstep,1) * Sim.gamma; %  + skew_delta_angle * ones(1,length(Sim.gamma)); % Skew_New
else % Skew_New
    Sim.skew_angle = theta_iniziale;   % delta theta per diff coenergia
end

%% Parametri simulazione
Sim.n = In.Vel;
Sim.sim_angle = In.Ang;

% Valori di corrente Id e Iq da simulare, cambiano se è presente lo skew
if (Sim.skew)
    Sim.Id_ = Sim.skew_Imod .* cos(Sim.skew_gamma*pi/180);
    Sim.Iq_ = Sim.skew_Imod .* sin(Sim.skew_gamma*pi/180);
else
    % Calcola anche Id e Iq per ogni caso
    Sim.Id_ = Sim.Id;
    Sim.Iq_ = Sim.Iq;
end

% 2013/06/10 MG modiciato il caricamento del file Parmachine, adesso è in RUN_SIM e si carica col nome del file.mn o della cartella 
%% Carica il caso ed imposta alcuni valori di default
% addpath ([cd SLASH 'M']);
% here = cd;
% 
% if not(exist('RQ','var'))
%     [FileName,FilePath] = uigetfile('*Machine.mat','Select the ParMachine Mat-file');
% end
% load([FilePath 'ParMachine.mat']);

% valori istantanei per ogni simulazione
Ft.legend = 't(ms) Fa(Vs) Fb(Vs) Fc(Vs) Ea(V) Eb(V) Ec(V) Id(Apk) Iq(Apk) Fd(Vs) Fq(Vs) T(Nm) CoEnergy(J)';
Ft.values = [];

MachineFileName = [FilePath Mac.MachineNameMn];

%% Elab
Sim.Nstep = In.Nstep;                             % Number of simulations (default: one step per airgap mesh step)
Sim.N_turn = Mac.N_turn;                          % spire in serie per fase
Sim.N_parallel = Mac.N_parallel;                  % numero di vie
Sim.N_cond = Mac.N_cond;                          % conduttori in Magnet

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Avvio delle simulazioni e salvataggi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MN6 = actxserver('Magnet.application');

RunS_RunCurve;

invoke(MN6,'close')

cd(here);
