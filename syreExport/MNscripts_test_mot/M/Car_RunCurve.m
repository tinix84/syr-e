%% VERSIONE 20 11 2011
% Car_RunCurve.m    - 20 06 08 - GMP

% called by AAB_RunIdentification.m
% given Id, Iq vectors runs "k x kk" MagNet simulations

% j: numeri complessi
clear j

% Sim: input di n x m simulazioni.
% CurveEstreme
%   - n = 4, m = numero punti per curva
% Matrice completa
%   - n = direzione secondaria, m = direzione principale

% Cas: input di 1 simulazione (caso.mn)
Cas.n          = Sim.n;
Cas.sim_angle  = Sim.sim_angle;
Cas.Nstep      = Sim.Nstep;
Cas.N_turn     = Sim.N_turn;
Cas.N_parallel = Sim.N_parallel;

%% Nome Cartella Backup
temp1 = num2str(Sim.Idmin,'%g');
temp2 = num2str(Sim.Idmax,'%g');
temp3 = num2str(Sim.Iqmin,'%g');
temp4 = num2str(Sim.Iqmax,'%g');

temp1(temp1 == '.') = 'A';
temp2(temp2 == '.') = 'A';
temp3(temp3 == '.') = 'A';
temp4(temp4 == '.') = 'A';

temp3=['Map_' temp1 '_' temp2 '_' temp3 '_' temp4];

if (In.CurveEstreme)
    temp3 = [temp3 '_4x_' num2str(In.Vel) 'rpm'];
else
    temp3 = [temp3 '_nxm_' num2str(In.Vel) 'rpm'];
end

%% Nome della Cartella completo di path
temp4 = [Mac.MachinePath 'CaseBackup' SLASH temp3]
temp5 = 'Y';

%% check SOVRASCRIVERE
% Y:        sovrascrive
% N:        break
% K:        continua vecchia mappa nella vecchia cartella

if isdir(temp4)
    temp5=questdlg('Caso già simulato. Sovrascrivere?','WARNING','Y','N','K','N');
end
if (temp5 == 'Y')
    mkdir([Mac.MachinePath 'CaseBackup' SLASH] , temp3);
elseif (temp5 == 'N')
    break
else % (Keep)
    warndlg('OK, completo la mappa','KEEP');
end

%% Ciclo for righe kk
for kk = 1:size(Sim.Id_,1)
    % Id Iq del caso (skew o no skew)
    Cas.Id_        = Sim.Id_(kk,:);
    Cas.Iq_        = Sim.Iq_(kk,:);

    Cas.Imod_      = abs(Cas.Id_ + j*Cas.Iq_);
    Cas.gamma_     = angle(Cas.Id_ + j*Cas.Iq_);

    %% Ciclo for colonne k
    for k = 1:size(Sim.Id_,2)
%         temporanea = [Cas.Id_(k) Cas.Iq_(k)]
        temp1 = num2str(Cas.Id_(k),'%g');
        temp2 = num2str(Cas.Iq_(k),'%g');

        temp1(temp1 == '.') = 'A';
        temp2(temp2 == '.') = 'A';

        %% Nome simulazione .mn (una simulazione per fetta skew)
        CaseFileName = [Mac.MachinePath 'CaseBackup' SLASH temp3 SLASH Mac.MachineName '_' temp1 '_' temp2 '.mn'];

        DocSV = invoke(MN6, 'openDocument',MachineFileName);

        % caso già girato (mappa interrotta o da ampliare)
%         return
        if (exist(CaseFileName,'file') == 2)
            save temp_Cas Cas
            load([CaseFileName(1:end-3) '.mat']);
            load temp_Cas
            % Append data mean (Fd Fq map)
            F_map.Fd(kk,k) = mean(Ft.values(:,10));
            F_map.Fq(kk,k) = mean(Ft.values(:,11));
            F_map.T(kk,k) = -mean(Ft.values(:,12));
%             keyboard 
        else
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

            % imposta in sequenza la simulazione eseguita
            Mac.caso = Mac.caso+1;
            Cas.caso = Mac.caso;

            % Run post process
            RunS_PostProcessing;

            % Append data mean (Fd Fq map)
            % valori medi
            F_map.Fd(kk,k) = Fluxd_0;
            F_map.Fq(kk,k) = Fluxq_0;
            F_map.T(kk,k) = Torque_0;

        end
        
        F_map.Id(kk,k) = Sim.Id_(kk,k);
        F_map.Iq(kk,k) = Sim.Iq_(kk,k);

        % Append data plots (all time values)
        Ft_map.values=[Ft_map.values;Ft.values];


        % save map (temp)
        save([temp4 '\map_F.mat'],'F_map');
        save([temp4 '\map_Ft.mat'],'Ft_map');
    end
end

%     % save map (final)
%     save([temp4 '\map_F.mat'],'F_map');
%     save([temp4 '\map_Ft.mat'],'Ft_map');
