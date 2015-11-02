%% VERSIONE 20 11 2011
% Caratterizzazione
% Input struct (In):
%   Vel - rpm
%   Ang -
%   skew
%   skew_angle
%   skew_Nstep
%   temperature
% Output:

switch_computer
load('temp_files\BATCH_SETUP.mat')

%% Parametri simulazione
Sim.n  = In.Vel;
Sim.sim_angle = In.Ang;
Sim.Idmin = In.Idmin;
Sim.Idmax = In.Idmax;
Sim.Iqmin = In.Iqmin;
Sim.Iqmax = In.Iqmax;
Sim.d_Id = In.d_Id;
Sim.d_Iq = In.d_Iq;

switch In.CurveEstreme
    case 0
        
        temp1 = In.Idmin : In.d_Id : In.Idmax;
        temp2 = In.Iqmin : In.d_Iq : In.Iqmax;
%         temp1=[0 10 30 60 100 120 140 160 180 200];
%         temp2=temp1;
        [Sim.Id_,Sim.Iq_] = meshgrid(temp1,temp2);
        
    case 1
        temp1 = In.Idmin : In.d_Id : In.Idmax;
        temp2 = In.Idmin : In.d_Id : In.Idmax;
        temp3 = In.Idmin * ones(size(temp1));
        temp4 = In.Idmax * ones(size(temp1));
        
        Sim.Id_ = [temp1;temp2;temp3;temp4];
        
        temp1 = In.Iqmin * ones(size(temp1));
        temp2 = In.Iqmax * ones(size(temp1));
        temp3 = In.Iqmin : In.d_Iq : In.Iqmax;
        temp4 = temp3;
        
        Sim.Iq_ = [temp1;temp2;temp3;temp4];
        
        clear temp1 temp2 temp3 temp4
        
    case 2
        
        temp1 = In.Idmin : In.d_Id : In.Idmax;
        Sim.Id_ = temp1;
        Sim.Iq_ = zeros(size(temp1));
        
    case 3
        
        temp1 = In.Iqmin : In.d_Iq : In.Iqmax;
        Sim.Id_ = zeros(size(temp1));
        Sim.Iq_ = temp1;
        
    otherwise
        
        return
        
end
% keyboard
%% Impostazioni per lo skewing
Sim.skew = In.skew;

if (Sim.skew)
    Sim.skew_Nstep = In.skew_Nstep;
    % Vettore skew_Imod
    %     Sim.skew_Imod = ones(Sim.skew_Nstep,1) * Sim.Imod;
    temp = In.skew_angle/Sim.skew_Nstep*(0:1:Sim.skew_Nstep-1);
    Sim.skew_angle = (temp - mean(temp))';   % Skew_New
    %     Sim.skew_gamma = ones(Sim.skew_Nstep,1) * Sim.gamma; %  + skew_delta_angle * ones(1,length(Sim.gamma)); % Skew_New
else
    Sim.skew_Nstep = 1;
    Sim.skew_angle = 0;   % Skew_New
end

% valori istantanei per ogni simulazione (salvato in RunS_PostProcessing.m)
Ft.legend = 't(ms) Fa(Vs) Fb(Vs) Fc(Vs) Ea(V) Eb(V) Ec(V) Id(Apk) Iq(Apk) Fd(Vs) Fq(Vs) T(Nm)';
Ft.values = [];

% valori istantanei, complessivo (salvato in Car_RunCurve.m)
Ft_map.legend = 't(ms) Fa(Vs) Fb(Vs) Fc(Vs) Ea(V) Eb(V) Ec(V) Id(Apk) Iq(Apk) Fd(Vs) Fq(Vs) T(Nm)';
Ft_map.values = [];

MachineFileName = [pathname MachineNameMn];

%% Elab
Sim.temperature = In.temperature;
Sim.Nstep = In.Nstep;     % Number of simulations (one step per airgap mesh step)
Sim.N_turn = Mac.N_turn;                          % spire in serie per fase
Sim.N_parallel = Mac.N_parallel;                  % numero di vie
Sim.N_cond = Mac.N_cond;                          % conduttori in Magnet


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Avvio delle simulazioni e salvataggi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Avvia le simulazioni
MN6 = actxserver('Magnet.application');
set(MN6, 'Visible', 1);

Car_RunCurve_SkewNew;

% invoke(MN6,'close')
Command='Call close(False)';
invoke(MN6, 'processCommand', Command);

cd(here)

