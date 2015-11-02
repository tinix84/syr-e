%% VERSIONE 09/06/2015
% incorpora RunSimulation e Caratterizzazione
% accetta vettori di velocità o temperatura (se uno è vettore l'altro deve
% essere scalare) e salva in cartelle diverse le diverse velocità
% (temperature)

clear all, close all, addpath M
answer = questdlg('Type of simulation','it''s time to choose','__SIMULATION__','IDENTIFICATION','__SIMULATION__');
load temp_files\last_for_default

switch answer
    case'__SIMULATION__'
               
        prompt={'SPEED RPM','delta position (elt deg)','n sim','skew yes or no', ...
            'skew angle (mech deg)','skew slices','temperature (°C)','Id vector','Iq vector (same size)','choose "SYRE" or "VDES" Input data format'};
        name='MOTOR SIMULATION SETUP';
        numlines=1;
        if length(setup) == 10
            defaultanswer={num2str(setup{1}),num2str(setup{2}),num2str(setup{3}),num2str(setup{4}),num2str(setup{5}), ...
                num2str(setup{6}),num2str(setup{7}),'50','50','SYRE'};
        else
            defaultanswer = setup;
        end
        setup=inputdlg(prompt,name,numlines,defaultanswer);
        
        temp_id = eval(setup{8}); % A
        temp_iq = eval(setup{9}); % A
        
        % Modulo di corrente in [A]  (puo' essere un vettore)
        In.Imod = abs(temp_id + 1i* temp_iq);
        
        % Fase della corrente in [gradi]  (puo' essere un vettore)
        In.gamma = angle(temp_id + 1i* temp_iq) * 180/pi;
        % In.gamma = 54.8;
        %% 03 06 09 - torque ripple - coenergia - test - (default = 0)
        theta_iniziale = 0;
        In.DataFormat=setup{10};
    case 'IDENTIFICATION'
        
        prompt={'SPEED RPM','delta position (elt deg)','n sim','skew yes or no', ...
            'skew angle (mech deg)','skew slices','temperature (°C)','Id min','Id max', ...
            'Iq min','Iq max','step Id','step Iq','Curve estreme','choose "SYRE" or "VDES" Input data format'};
        name='MOTOR IDENTIFICATION SETUP';
        numlines=1;

        if length(setup) == 15
            defaultanswer = setup;
        else
            defaultanswer={num2str(setup{1}),num2str(setup{2}),num2str(setup{3}),num2str(setup{4}),num2str(setup{5}), ...
                num2str(setup{6}),num2str(setup{7}),'0','100','0','100','10','10','1','SYRE'};
        end
        
        setup=inputdlg(prompt,name,numlines,defaultanswer);
        
        % estremi di corrente
        In.Idmin        =   eval(setup{8});          % A
        In.Idmax        =	eval(setup{9});        % A
        In.Iqmin        =	eval(setup{10});          % A
        In.Iqmax        =	eval(setup{11});        % A
        % passo
        In.d_Id         =	eval(setup{12});        % A
        In.d_Iq         =	eval(setup{13});         % A
        % (0) Matrice mxn, (1) CurveEstreme 4x, (2) singola curva lungo id , (3) singola curva lungo iq
        In.CurveEstreme = eval(setup{14});
        In.DataFormat=setup{15};
    otherwise
        disp('Unknown method.')
end

%% DATI GENERALI
% Velocita di rotazione del motore in [rpm]
In.Vel = eval(setup{1});
% Periodo di simulazione della macchina in [elettrici]
In.Ang  = eval(setup{2});
% Numero di passi (default 1 passo = 0.5 gradi elettrici)
In.Nstep = eval(setup{3});
% Skewing (1->si) (0->no) skewing con i dati di angolo e fette
In.skew = eval(setup{4});
% Angolo di skew in GRADI MECCANICI % Skew_New
In.skew_angle = eval(setup{5});
% Numero di fette
In.skew_Nstep = eval(setup{6});
% temperatura (°C)
In.temperature = eval(setup{7});
%%
    prompt={'Stator Iron used', 'rotor iron used','Magnet'};
    name='IRON and MAGNET IDENTIFICATION';
    numlines=1;
    defaultanswer={'M250-35A','M250-35A','None'};
    IronName=inputdlg(prompt,name,numlines,defaultanswer);
%% Variabili da modificare in Batch
if length(In.temperature) > 1
    RQ.name = 'temperature';
    RQ.values = In.temperature;
end

if length(In.Vel) > 1
    RQ.name = 'Vel';
    RQ.values = In.Vel;
end

here = cd;
[FileName,FilePath] = uigetfile([FilePath '**.mat'],'Select the machine dat file Mat-file');
load([FilePath,FileName]);

if (strcmp(In.DataFormat,'SYRE'))
   Mac=genMacDataFromGeo(geo);
   Mac.n_mag_simulati=0;
end

% machine name for load and saving
   MachineNameMn=[FileName(1:end-4),'.mn'];
   MachineName=FileName(1:end-4); 

if not(exist('RQ','var'))
    RQ.name = 'temperature';
    RQ.values = In.temperature;
end

save('temp_files\last_for_default','setup','FilePath');
%     keyboard
for batch_step = 1:length(RQ.values)
    
%     warndlg(['running batch - step' num2str(batch_step) 'of ' num2str(length(RQ.values))]);
    eval(['In.' RQ.name ' = ' num2str(RQ.values(batch_step)) ';']);
    save('temp_files\BATCH_SETUP.mat','In','IronName')
    
    if isfield(In,'CurveEstreme')
        AAB_Caratterizzazione_SkewNew_with_batch
    else
        AAA_RunSimulation_SkewNew_with_batch
    end
    
end



