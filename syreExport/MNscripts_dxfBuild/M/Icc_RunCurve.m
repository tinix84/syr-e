%% VERSIONE 20 11 2011
% Icc_RunCurve.m    - 19 12 08 - GMP

% called by AAC_Find_Icc.m
% runs iterative MagNet simulations to find idq = 0 + j*icc

% j: numeri complessi

clear j
% Cas: vettore del singolo caso     (1 x 1)
% Sim: append n simulazioni         (1 x n)

Sim.n          = Cas.n;
Sim.N_turn     = Cas.N_turn;
Sim.N_parallel = Cas.N_parallel;
Sim.N_cond       = Cas.N_cond;
Sim.temperature = Cas.temperature;

% nome cartella di backup
temp3 = 'Icc';

%% Nome della Cartella completo di path
temp4 = [Mac.MachinePath 'CaseBackup' SLASH temp3];
temp5 = 'Y';

%% check SOVRASCRIVERE
% Y:        sovrascrive
% N:        break
% K:        continua nella vecchia cartella

if isdir(temp4)
    temp5=questdlg('Caso già simulato. Sovrascrivere?','WARNING','Y','N','K','N');
end
if (temp5 == 'Y')
    mkdir([Mac.MachinePath 'CaseBackup' SLASH] , temp3);
elseif (temp5 == 'N')
    break
else % (Keep)
    warndlg('OK, completo il lavoro','KEEP');
    temp3 = [temp3 '_refine'];
    temp4 = [Mac.MachinePath 'CaseBackup' SLASH temp3];
    mkdir([Mac.MachinePath 'CaseBackup' SLASH] , temp3);
end

k = 1;
% iq_min = 0;
% iq_max = 100;
Fluxq_0 = 1;

%% Ciclo for kk
while (abs(Fluxq_0) > 0.001) && (k < 10)
    
    % Cas (cas_xxx.mn)
    Sim.sim_angle(k)  = Cas.sim_angle;
    Sim.Nstep(k)      = Cas.Nstep;

    Cas.Imod_(k)   = abs(Cas.Id + j*Cas.Iq);
    Cas.gamma_(k)  = angle(Cas.Id + j*Cas.Iq);
    Cas.Id_(k)     = Cas.Id;
    Cas.Iq_(k)     = Cas.Iq;


    %% Ciclo for colonne k
    %     for k = 1:size(Sim.Id_,2)
    %         temporanea = [Cas.Id_(k) Cas.Iq_(k)]
    temp1 = num2str(Cas.Id,'%g');
    temp2 = num2str(Cas.Iq,'%g');

    temp1(temp1 == '.') = 'A';
    temp2(temp2 == '.') = 'A';

    %% Nome simulazione .mn (una simulazione per fetta skew)
    CaseFileName = [Mac.MachinePath 'CaseBackup' SLASH temp3 SLASH Mac.MachineName '_' temp1 '_' temp2 '.mn'];

    DocSV = invoke(MN6, 'openDocument',MachineFileName);

    % caso già girato (mappa interrotta o da ampliare)
    %     if (exist(CaseFileName,'file') == 2)
    %         save temp_Cas Cas
    %         load([CaseFileName(1:end-3) '.mat']);
    %         load temp_Cas
    %         % Append data mean (Fd Fq map)
    %         F_icc.Fd(kk,k) = mean(Ft.values(:,10));
    %         F_icc.Fq(kk,k) = mean(Ft.values(:,11));
    %         F_icc.T(kk,k) = -mean(Ft.values(:,12));
    %         %             keyboard
    %     else
    %% caso nuovo, lancia simulazione
    Doc = invoke(MN6, 'getDocument');
    View = invoke(Doc, 'getCurrentView');
    Solution = invoke(Doc, 'getSolution');
    calculator = invoke (MN6, 'getStackCalculator');

    % set the operating parameters
    RunS_SetParCase;

    % solve
    invoke(Doc, 'solveTransient2dWithMotion');

    %save the model (.mn)
    DocSV = invoke(MN6, 'saveDocument',CaseFileName);

    %         DocSV = invoke(MN6, 'saveDocument','testa.mn');
    % imposta in sequenza la simulazione eseguita
    %         Mac.caso = Mac.caso+1;
    %         Cas.caso = Mac.caso;

    % Run post process
    RunS_PostProcessing;

    % Append data mean (Fd Fq map)
    % valori medi
    F_icc.Fd(k) = Fluxd_0;
    F_icc.Fq(k) = Fluxq_0;
    F_icc.T(k) = Torque_0;

    %     end

    F_icc.Id(k) = Cas.Id;
    F_icc.Iq(k) = Cas.Iq;

    % Append data plots (all time values)
    Ft_icc.values=[Ft_icc.values;Ft.values];


    % save map (temp)
    save([temp4 '\Icc_F.mat'],'F_icc');
    save([temp4 '\Icc_Ft.mat'],'Ft_icc');
    %     end

    Sim.Id(k) = Cas.Id;
    Sim.Iq(k) = Cas.Iq;

    if (Fluxq_0 > 0)
        % ridurre
        iq_max = Cas.Iq;
    else
        % aumentare
        iq_min = Cas.Iq;
    end

    iq_min
    Cas.Iq = (iq_min + iq_max)/2
    iq_max
    
    k = k +1;
end

