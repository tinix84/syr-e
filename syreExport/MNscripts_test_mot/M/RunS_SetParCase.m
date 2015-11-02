%% VERSIONE 20 11 2011
% RunS_SetParCase.m
% ..
% ..

%% type of solution (order)
% Polynomial Order = 1 (3 field points per triangle)
% Command = ['Call getDocument().setParameter("", "PolynomialOrder", "1", infoNumberParameter)'];
% invoke(MN6, 'processCommand', Command);

% SourcesOnAtTransientStart: smooth field at t = 0
Command = ['Call getDocument().setParameter("", "SourcesOnAtTransientStart", "Yes", infoStringParameter)'];
invoke(MN6, 'processCommand', Command);

% 09 03 2010 - MeshAllowUnconstrainedHoles
Command = ['Call getDocument().setParameter("", "MeshAllowUnconstrainedHoles", "Yes", infoStringParameter)'];
invoke(MN6, 'processCommand', Command);

% Model temperature
Command = ['Call getDocument().setParameter("", "Temperature", "' num2str(Sim.temperature) '", infoNumberParameter)'];
invoke(MN6, 'processCommand', Command);

% set transient options (time0 etc..)
time_end = 1000 * Cas.sim_angle/(360*Mac.p)*60/Cas.n;   % ms
invoke(Doc,'setFixedIntervalTimeSteps',0, time_end, time_end/Cas.Nstep);

% Iron loss time interval - global
Command = ['Call getDocument().setParameter("", "TransientAveragePowerLossStartTime", "0", infoNumberParameter)'];
invoke(MN6, 'processCommand', Command);

Command = ['Call getDocument().setParameter("", "TransientAveragePowerLossEndTime", "' num2str(time_end) ' %ms", infoNumberParameter)'];
invoke(MN6, 'processCommand', Command);

% Iron loss time interval - rotor
% time_start_pferot = time_end - time_end/Mac.Q;
%% temp -- macchina vagati
% time_start_pferot = time_end - Mac.q/Mac.Q*time_end;
%% 2013/11/11 MG valutazione delle perdite rotoriche sugli ultimi 60°
% time_start_pferot = time_end - 60/(360*Mac.p)*60/Cas.n;
% 
% Command = ['Call getDocument().setParameter("rotor", "TransientAveragePowerLossStartTime", "' num2str(time_start_pferot) ' %ms", infoNumberParameter)'];
% invoke(MN6, 'processCommand', Command);
% 
% Command = ['Call getDocument().setParameter("rotor", "TransientAveragePowerLossEndTime", "' num2str(time_end) ' %ms", infoNumberParameter)'];
% invoke(MN6, 'processCommand', Command);
%%  %%%%%%%%%%
Command = ['Call getDocument().setMotionSourceType("Motion#1", infoVelocityDriven)'];
invoke(MN6, 'processCommand', Command);
%%

% Motion#1: transient angle Vs time
% start time = 0
Command = 'REDIM ArrayOfValues1(1)';
invoke(MN6, 'processCommand', Command);
Command = 'ArrayOfValues1(0)= 0';
invoke(MN6, 'processCommand', Command);
% end time == sim_angle (elt deg)
Command = ['ArrayOfValues1(1)= ' num2str(time_end)];
invoke(MN6, 'processCommand', Command);

% start angle = Cas.skew_angle
Command = 'REDIM ArrayOfValues2(1)';
invoke(MN6, 'processCommand', Command);
Command = ['ArrayOfValues2(0)= ' num2str(Cas.skew_angle)]; % Skew_New
invoke(MN6, 'processCommand', Command);
% end angle = start + sim_angle
Command = ['ArrayOfValues2(1)= ' num2str(Cas.sim_angle/Mac.p + Cas.skew_angle)];   %sim_angle elettrici % Skew_New
invoke(MN6, 'processCommand', Command);
Command = 'Call getDocument().setMotionPositionVsTime("Motion#1", ArrayOfValues1, ArrayOfValues2)';
invoke(MN6, 'processCommand', Command);

%%
phase_index = {'U','V','W'};
if isfield(Mac,'th0')
    % casi successivi feb 2010 (+90 perchè Waveform in Magnet\coil è un sin)
    phase_angle = [0 -120 120] + (Cas.gamma_(k)*180/pi + Mac.th0 + 90);
else
    % casi precedenti feb 2010
    phase_angle = [0 -120 120] - 180/pi*(pi-Cas.gamma_(k));
end
for ii = 1:3
    Command=['Call getDocument().setCoilSourceType("',phase_index{ii},'", infoCurrentDriven)'];
    invoke(MN6, 'processCommand', Command);

    Command = 'REDIM ArrayOfValues(5)';
    invoke(MN6, 'processCommand', Command);
    Command = 'ArrayOfValues(0)= 0';
    invoke(MN6, 'processCommand', Command);
    Command = ['ArrayOfValues(1)= ' num2str(Cas.Imod_(k)/Cas.N_parallel)];
    invoke(MN6, 'processCommand', Command);
    Command = ['ArrayOfValues(2)= ' num2str(Cas.n/60*Mac.p)];
    invoke(MN6, 'processCommand', Command);
    Command = 'ArrayOfValues(3)= 0';
    invoke(MN6, 'processCommand', Command);
    Command = 'ArrayOfValues(4)= 0';
    invoke(MN6, 'processCommand', Command);
    Command = ['ArrayOfValues(5)= ' num2str(phase_angle(ii))];
    invoke(MN6, 'processCommand', Command);
    Command = ['Call getDocument.setSourceWaveform("' phase_index{ii} '","SIN", ArrayOfValues)'];
    invoke(MN6, 'processCommand', Command);
    % Command = ['getDocument.setParameter("' index{i} '", "WaveFormValues", "[0, ' num2str(Imod) ', ' num2str(n/60*p) ', 0, 0, 90]", infoArrayParameter)'];
    % invoke(MN6, 'processCommand', Command);
end

%numbers of turns
Command = ['Call getDocument().setCoilNumberOfTurns("U", ' num2str(Sim.N_cond) ')'];
invoke(MN6, 'processCommand', Command);
Command = ['Call getDocument().setCoilNumberOfTurns("V", ' num2str(Sim.N_cond) ')'];
invoke(MN6, 'processCommand', Command);
Command = ['Call getDocument().setCoilNumberOfTurns("W", ' num2str(Sim.N_cond) ')'];
invoke(MN6, 'processCommand', Command);

