% function [SOL] = simulate_xdeg_MN(geo,per,mat,eval_type,pathname,filename)
function [SOL] = simulate_xdeg_MN(geo,per,eval_type,pathname,filename,h)

% number of simulation that must be done respect to eval type
switch eval_type
    %     case 'MO_OA'
    %         gamma = gamma_in;
    %         nsim = geo.nsim_MOOA;
    %         xdeg = geo.delta_sim_MOOA;
    %         randFactor = geo.randFactor;
    %     case 'MO_GA'
    %         gamma = gamma_in;
    %         nsim = geo.nsim_MOOA;
    %         xdeg = geo.delta_sim_MOOA;
    %         randFactor = geo.randFactor;
    case 'singt'
        nsim = geo.nsim_singt;
        xdeg = geo.delta_sim_singt;
        gamma = per.gamma;
    case 'singm'
        nsim = geo.nsim_singt;
        xdeg = geo.delta_sim_singt;
        gamma = per.gamma;
    otherwise
        error('MagNet Simulations not available during optimization!');
end

randFactor = 0;


% pathname = pwd();

th0 = geo.th0;
p   = geo.p;
r  = geo.r;
gap = geo.g;
ns  = geo.ns;
pc  = 360/(ns*p)/2;
ps  = geo.ps;
n3phase = geo.n3phase; %AS number of 3-phase circuits
Q = geo.ns*geo.p;
skew_angle = 0;
N_parallel = 1;
N_cond = geo.Ns/geo.p/geo.q/size(geo.avv,1);
Br = per.BrPP;
tempPP = per.tempPP;
q = geo.q;
gamma_ = gamma*pi/180;

% simulation angle
gradi_da_sim=180/p*ps;

% Hc = per.BrPP/(4e-7*pi);

%SOL = [];

% evaluation of the phase current values for all positions to be simulated
iAmp = per.overload*calc_io(geo,per);
% iAmpCoil = iAmp*geo.Nbob*geo.n3phase; %AS

id = iAmp * cos(gamma * pi/180);
iq = iAmp * sin(gamma * pi/180);



%load('temp_files\BATCH_SETUP.mat')

Imod = abs(id + 1i* iq);       % Modulo corrente (A)

%% Avvio delle simulazioni e salvataggi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RunS_RunCurve

% clear j


%     %% Nome Cartella (1 Cartella per ogni Caso)  !!!!cambiare!!!!!!
%     temp1 = num2str(Imod,'%0.2f');
%     temp1(temp1 == '.') = 'A';
%     temp2 = num2str(gamma,'%0.1f');
%     temp2(temp2 == '.') = 'G';
%     temp3=[temp1 '_' temp2 '_' num2str(xdeg) 'degs'];
%
% %     %% Cartella del Caso (Nstep simulazioni)
%
%     SLASH='\';
%     temp4 = [pathname,FileName(1:end-4) 'CaseBackup' SLASH temp3 '_' num2str(xdeg) 'rpm'];
%
%     if exist('RQ','var')
%         %% batch
%         temp4 = [temp4 '_' RQ.name '_' num2str(RQ.values(batch_step))];
%     end
%     temp5 = 'Yes';
%     if isdir(temp4)
%         temp5=questdlg('Caso già simulato. Sovrascrivere?','WARNING');
%     end
%     if (temp5 == 'Yes')
%         mkdir(temp4);
%     end
%     CaseDir = temp4;
%     clear temp4 temp5


if  (eval_type == 'singt')
    iStr=num2str(Imod,3); iStr = strrep(iStr,'.','A');
    gammaStr=num2str(gamma,4); gammaStr = strrep(gammaStr,'.','d');
    if isempty(strfind(gammaStr, 'd'))
        gammaStr = [gammaStr 'd'];
    end
    
    CaseDir = [pathname, filename(1:end-4) '_T_eval_',iStr,'_',gammaStr '_MN'];
    mkdir(CaseDir);
else
    Idstr=num2str(id,3); Idstr = strrep(Idstr,'.','A');
    Iqstr=num2str(iq,3); Iqstr = strrep(Iqstr,'.','A');
    
    if isempty(strfind(Idstr, 'A'))
        Idstr = [Idstr 'A'];
    end
    if isempty(strfind(Iqstr, 'A'))
        Iqstr = [Iqstr 'A'];
    end
    
    CaseDir = [pathname 'CaseBackup\', filename(1:end-4),'_' Idstr 'x' Iqstr '_MN'];
    mkdir(CaseDir );
    
end
%%

CaseFileName = [CaseDir '\' filename(1:end-4) '.mn'];  % folder save .mn
filename = [filename(1:end-4) '.mn'];

% MN6 = actxserver('Magnet.application');
% set(MN6, 'Visible', 1);
% DocSV = invoke(MN6, 'openDocument',[pathname filename]);
% Doc = invoke(MN6, 'getDocument');
% View = invoke(Doc, 'getCurrentView');
% Solution = invoke(Doc, 'getSolution');
% calculator = invoke (MN6, 'getStackCalculator');
% h = OpenMagnet(1);
DocSV = invoke(h.magnetHandler, 'openDocument',[pathname filename]);
Doc = invoke(h.magnetHandler, 'getDocument');
View = invoke(Doc, 'getCurrentView');
Solution = invoke(Doc, 'getSolution');
calculator = invoke (h.magnetHandler, 'getStackCalculator');

% set the operating parameters
if (Br==0)
    n_mag_simulati=0;
else
    %n_mag_simulati=size(geo.BLKLABELS.rotore.BarName,1);
    tmp=length(geo.BLKLABELS.rotore.xy(3,:));
    n_mag_simulati=length(tmp(tmp==6));
end
%%
% RunS_SetParCase.m
% type of solution (order)
% Polynomial Order = 1 (3 field points per triangle)
% Command = ['Call getDocument().setParameter("", "PolynomialOrder", "1", infoNumberParameter)'];
% invoke(MN6, 'processCommand', Command);

% SourcesOnAtTransientStart: smooth field at t = 0
Command = ['Call getDocument().setParameter("", "SourcesOnAtTransientStart", "Yes", infoStringParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

% 09 03 2010 - MeshAllowUnconstrainedHoles
Command = ['Call getDocument().setParameter("", "MeshAllowUnconstrainedHoles", "Yes", infoStringParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

% Model temperature
Command = ['Call getDocument().setParameter("", "Temperature", "' num2str(tempPP) '", infoNumberParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

% set transient options (time0 etc..)
time_end = 1000 * xdeg/(360*p)*60/per.EvalSpeed;   % ms
invoke(Doc,'setFixedIntervalTimeSteps',0, time_end, time_end/nsim);

% Iron loss time interval - global
Command = ['Call getDocument().setParameter("", "TransientAveragePowerLossStartTime", "0", infoNumberParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

Command = ['Call getDocument().setParameter("", "TransientAveragePowerLossEndTime", "' num2str(time_end) ' %ms", infoNumberParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

% Iron loss time interval - rotor
time_start_pferot = time_end - time_end/Q;
%% temp -- macchina vagati
time_start_pferot = time_end - q/Q*time_end;
Command = ['Call getDocument().setParameter("rotor", "TransientAveragePowerLossStartTime", "' num2str(time_start_pferot) ' %ms", infoNumberParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

Command = ['Call getDocument().setParameter("rotor", "TransientAveragePowerLossEndTime", "' num2str(time_end) ' %ms", infoNumberParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);
%% MATTEO  no trovato il comando per velocity driven inerito in questa  parte VERIFICARE:
Command = ['Call getDocument().setMotionSourceType("Motion#1", infoVelocityDriven)'];
invoke(h.magnetHandler, 'processCommand', Command);
%%

% Motion#1: transient angle Vs time
% start time = 0
Command = 'REDIM ArrayOfValues1(1)';
invoke(h.magnetHandler, 'processCommand', Command);
Command = 'ArrayOfValues1(0)= 0';
invoke(h.magnetHandler, 'processCommand', Command);
% end time == sim_angle (elt deg)
Command = ['ArrayOfValues1(1)= ' num2str(time_end)];
invoke(h.magnetHandler, 'processCommand', Command);

% start angle = Cas.skew_angle
Command = 'REDIM ArrayOfValues2(1)';
invoke(h.magnetHandler, 'processCommand', Command);
Command = ['ArrayOfValues2(0)= ' num2str(skew_angle)]; % Skew_New
invoke(h.magnetHandler, 'processCommand', Command);
% end angle = start + sim_angle
Command = ['ArrayOfValues2(1)= ' num2str(xdeg/p + skew_angle)];   %sim_angle elettrici % Skew_New
invoke(h.magnetHandler, 'processCommand', Command);
Command = 'Call getDocument().setMotionPositionVsTime("Motion#1", ArrayOfValues1, ArrayOfValues2)';
invoke(h.magnetHandler, 'processCommand', Command);

%%
phase_index = {'U','V','W'};
% if isfield(Mac,'th0')
if exist('th0')
    % casi successivi feb 2010 (+90 perchè Waveform in Magnet\coil è un sin)
    phase_angle = [0 -120 120] + (gamma_*180/pi + th0 + 90);
else
    % casi precedenti feb 2010
    phase_angle = [0 -120 120] - 180/pi*(pi-gamma_);
end
for ii = 1:3
    Command=['Call getDocument().setCoilSourceType("',phase_index{ii},'", infoCurrentDriven)'];
    invoke(h.magnetHandler, 'processCommand', Command);
    
    Command = 'REDIM ArrayOfValues(5)';
    invoke(h.magnetHandler, 'processCommand', Command);
    Command = 'ArrayOfValues(0)= 0';
    invoke(h.magnetHandler, 'processCommand', Command);
    Command = ['ArrayOfValues(1)= ' num2str(Imod/N_parallel)];
    invoke(h.magnetHandler, 'processCommand', Command);
    Command = ['ArrayOfValues(2)= ' num2str(per.EvalSpeed/60*p)];
    invoke(h.magnetHandler, 'processCommand', Command);
    Command = 'ArrayOfValues(3)= 0';
    invoke(h.magnetHandler, 'processCommand', Command);
    Command = 'ArrayOfValues(4)= 0';
    invoke(h.magnetHandler, 'processCommand', Command);
    Command = ['ArrayOfValues(5)= ' num2str(phase_angle(ii))];
    invoke(h.magnetHandler, 'processCommand', Command);
    Command = ['Call getDocument.setSourceWaveform("' phase_index{ii} '","SIN", ArrayOfValues)'];
    invoke(h.magnetHandler, 'processCommand', Command);
    % Command = ['getDocument.setParameter("' index{i} '", "WaveFormValues", "[0, ' num2str(Imod) ', ' num2str(n/60*p) ', 0, 0, 90]", infoArrayParameter)'];
    % invoke(h.magnetHandler, 'processCommand', Command);
end

%numbers of turns
Command = ['Call getDocument().setCoilNumberOfTurns("U", ' num2str(N_cond) ')'];
invoke(h.magnetHandler, 'processCommand', Command);
Command = ['Call getDocument().setCoilNumberOfTurns("V", ' num2str(N_cond) ')'];
invoke(h.magnetHandler, 'processCommand', Command);
Command = ['Call getDocument().setCoilNumberOfTurns("W", ' num2str(N_cond) ')'];
invoke(h.magnetHandler, 'processCommand', Command);
% keyboard
%%
% solve
invoke(Doc, 'solveTransient2dWithMotion');

RunS_PostProcessingMN;
%   RunS_PostProcessing;

%save the model (.mn)
Command=['Call getDocument().save("',CaseFileName,'")'];
invoke(h.magnetHandler, 'processCommand', Command);


%% output
tempo   = time(2:end);         %ms
th = (geo.th0 * pi/180 + tempo/1000 * per.EvalSpeed * pi/30 * geo.p);

SOL.pm_loss = pm_loss;
SOL.Ppm = Ppm;
if (xdeg == 360)
    SOL.Pfes_h = Pfes_h;
    SOL.Pfes_c = Pfes_c;
    SOL.Pfer_h = Pfer_h;
    SOL.Pfer_c = Pfer_c;
    SOL.PjrBar = PjrBar;
end
SOL.th = (th)';
SOL.id = tmp1(2:end)';
SOL.iq = tmp2(2:end)';
SOL.fd = Fluxd';
SOL.fq = Fluxq';
SOL.T = (mean([-Torque_1sim(:,1) Torque_1sim(:,2)],2))';

% Command='Call close(False)';
% invoke(h.magnetHandler, 'processCommand', Command);

cd(cd);






