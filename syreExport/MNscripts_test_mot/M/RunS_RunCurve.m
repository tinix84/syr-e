%% VERSIONE 20 11 2011
% RunCurve
% given Id, Iq vectors runs "k x kk" MagNet simulations

% j: numeri complessi
clear j

% Sim: input di n x m simulazioni. n = Nstep skew, m = numero di punti di lavoro
% Cas: input di 1 simulazione (caso.mn)
Cas.n          = Sim.n;
Cas.sim_angle  = Sim.sim_angle;
Cas.Nstep      = Sim.Nstep;
Cas.N_turn     = Sim.N_turn;
Cas.N_parallel = Sim.N_parallel;
Cas.skew       = Sim.skew;
Cas.skew_Nstep = Sim.skew_Nstep;

%% Ciclo for Casi
for kk = 1:size(Sim.Id_,2)
    % Id Iq del caso (skew o no skew)
    Cas.Imod       = Sim.Imod(kk);
    Cas.gamma      = Sim.gamma(kk);
    Cas.Id_        = Sim.Id_(:,kk);
    Cas.Iq_        = Sim.Iq_(:,kk);
    if (Sim.skew)
        Cas.Imod_      = Sim.skew_Imod(:,kk);
        Cas.gamma_     = Sim.skew_gamma(:,kk)*pi/180;
    else
        Cas.Imod_      = Cas.Imod;
        Cas.gamma_     = Cas.gamma*pi/180;
    end
    %% Nome Cartella (1 Cartella per ogni Caso)
    temp1 = num2str(Cas.Imod,'%0.2f');
    temp1(temp1 == '.') = 'A';
    temp2 = num2str(Cas.gamma,'%0.1f');
    temp2(temp2 == '.') = 'G';
    temp3=[temp1 '_' temp2 '_' num2str(Cas.sim_angle) 'degs'];
    if (Cas.skew)
        temp3 = [temp3 '_skew'];
    end
    
    %% Cartella del Caso (Nstep simulazioni)
    %     temp4 = [FilePath 'CaseBackup' SLASH temp3];
    % 2013/06/10 MG modificata la cartella di salvataggio dei dati,
    % adesso il nome è lo stesso del file ParMachine usato.
    temp4 = [pathname,FileName(1:end-4) 'CaseBackup' SLASH temp3 '_' num2str(Cas.n) 'rpm'];
    
    if exist('RQ','var')
        %% batch
        temp4 = [temp4 '_' RQ.name '_' num2str(RQ.values(batch_step))];
    end
    temp5 = 'Yes';
    if isdir(temp4)
        temp5=questdlg('Caso già simulato. Sovrascrivere?','WARNING');
    end
    if (temp5 == 'Yes')
        mkdir(temp4);
    else
        break
    end
    CaseDir = temp4;
    clear temp4 temp5
    %% Ciclo for: fette skew
    for k = 1:size(Cas.Id_,1)
        Cas.skew_angle = Sim.skew_angle(k);  % Skew_New
        temp4 = ['_' num2str(k)];
        %% Nome simulazione .mn (una simulazione per fetta skew)
        CaseFileName = [CaseDir SLASH MachineName temp4 '.mn'];
        
        set(MN6, 'Visible', 1);
        DocSV = invoke(MN6, 'openDocument',MachineFileName);
        Doc = invoke(MN6, 'getDocument');
        View = invoke(Doc, 'getCurrentView');
        Solution = invoke(Doc, 'getSolution');
        calculator = invoke (MN6, 'getStackCalculator');
        
        % set the operating parameters
        RunS_SetParCase;
        
        % solve
        invoke(Doc, 'solveTransient2dWithMotion');
        % imposta in sequenza la simulazione eseguita
        Mac.caso = Mac.caso+1;
        Cas.caso = Mac.caso;
        
        % Run post process
        RunS_PostProcessing;
        
        %save the model (.mn)
        Command=['Call getDocument().save("',CaseFileName,'")'];
        invoke(MN6, 'processCommand', Command);
        
    end
end
clear temp1 temp2 temp3 temp4