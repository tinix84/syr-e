%% VERSIONE 21 11 2011
% Car_RunCurve_SkewNew.m    - 14 04 09 - GMP

% called by AAB_Caratterizzazione_SkewNew.m
% - given Id, Iq vectors runs "k x kk x 1" MagNet simulations
% - if (Sim.skew), runs "k x kk x kkk" simulations

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

switch In.CurveEstreme
    case 0
        temp3 = [temp3 '_nxm_' num2str(In.Vel) 'rpm'];
    case 1
        temp3 = [temp3 '_4x_' num2str(In.Vel) 'rpm'];
    case 2
        temp3 = [temp3 '_d1x_' num2str(In.Vel) 'rpm'];
    case 3
        temp3 = [temp3 '_q1x_' num2str(In.Vel) 'rpm'];
    otherwise
        return
end

%% Nome della Cartella completo di path
temp4 = [pathname,FileName(1:end-14) 'CaseBackup' SLASH temp3];
% if exist('BATCH_SETUP.mat')
if exist('RQ','var')
    temp4 = [temp4 '_' RQ.name '_' num2str(RQ.values(batch_step))];
end
temp5 = 'Y';

%% check SOVRASCRIVERE
% Y:        sovrascrive
% N:        break
% K:        continua vecchia mappa nella vecchia cartella
if isdir(temp4)
    temp5=questdlg('Caso già simulato. Sovrascrivere?','WARNING','Y','N','K','N');
end
if (temp5 == 'Y')
    mkdir(temp4)
elseif (temp5 == 'N')
    return
else % (Keep)
    warndlg('OK, completo la mappa','KEEP');
end
% keyboard
%% Ciclo for righe kk
for kk = 1:size(Sim.Id_,1)
    % Id Iq del caso (skew o no skew)
    Cas.Id_        = Sim.Id_(kk,:);
    Cas.Iq_        = Sim.Iq_(kk,:);
    
    Cas.Imod_      = abs(Cas.Id_ + 1j*Cas.Iq_);
    Cas.gamma_     = angle(Cas.Id_ + 1j*Cas.Iq_);
    
    %% Ciclo for colonne k
    for k = 1:size(Sim.Id_,2)
        %         temporanea = [Cas.Id_(k) Cas.Iq_(k)]
        temp1 = num2str(Cas.Id_(k),'%g');
        temp2 = num2str(Cas.Iq_(k),'%g');
        
        temp1(temp1 == '.') = 'A';
        temp2(temp2 == '.') = 'A';
        
        for kkk = 1:Sim.skew_Nstep
            Cas.skew_angle  = Sim.skew_angle(kkk);
            %% Nome simulazione .mn (una simulazione per fetta skew)
            if (kkk == 1)
                CaseFileName = [temp4 SLASH MachineName '_' temp1 '_' temp2 '.mn'];
            else
                temp6 = num2str(kkk);
                CaseFileName = [temp4 SLASH MachineName '_' temp1 '_' temp2 '_' temp6 '.mn'];
            end
            
            DocSV = invoke(MN6, 'openDocument',MachineFileName);
            
            % caso già girato (mappa interrotta o da ampliare)
            if (exist(CaseFileName,'file') == 2)
                save temp_Cas Cas Sim
                load([CaseFileName(1:end-3) '.mat']);
                load temp_Cas
            else
                %% caso nuovo, lancia simulazione
                Doc = invoke(MN6, 'getDocument');
                View = invoke(Doc, 'getCurrentView');
                Solution = invoke(Doc, 'getSolution');
                calculator = invoke (MN6, 'getStackCalculator');
%                 keyboard
                % set the operating parameters
                RunS_SetParCase;
                
                % solve
                invoke(Doc, 'solveTransient2dWithMotion');
                
                % imposta in sequenza la simulazione eseguita
                Mac.caso = Mac.caso+1;
                Cas.caso = Mac.caso;
                
                % Run post process
                RunS_PostProcessing;

            end
            % Append data mean (Fd Fq map)
            F_map.Fd(kk,k,kkk) = mean(Ft.values(:,10));
            F_map.Fq(kk,k,kkk) = mean(Ft.values(:,11));
            F_map.T(kk,k,kkk) = -mean(Ft.values(:,12));
            F_map.dT(kk,k,kkk) = std(Ft.values(:,12));
            if Sim.sim_angle == 360
                F_map.Pfes_h(kk,k,kkk) = Pfes_h;
                F_map.Pfes_c(kk,k,kkk) = Pfes_c;
                F_map.Pfer_h(kk,k,kkk) = Pfer_h;
                F_map.Pfer_c(kk,k,kkk) = Pfer_c;
                F_map.Ppm(kk,k,kkk) = Ppm;
                F_map.PjrBar(kk,k,kkk)=PjrBar;
                F_map.stator_iron=IronName(1);
                F_map.rotor_iron=IronName(2);
                F_map.PMname=IronName(3);
                F_map.speed=In.Vel;
            end
            F_map.Id(kk,k,kkk) = Sim.Id_(kk,k);
            F_map.Iq(kk,k,kkk) = Sim.Iq_(kk,k);
            
        end     % skew
        
        % Append data plots (all time values)
        Ft_map.values=[Ft_map.values;Ft.values];
        
        % save map (temp)
        save([temp4 '\map_F.mat'],'F_map');
%         save([temp4 '\map_Ft.mat'],'Ft_map');
%         save([temp4 '\Fmap_Loss.mat'],'Fmap_Loss','stator_iron','rotor_iron','PMname');

    end     % colonne
end         % righe

% DocSV = invoke(MN6, 'saveDocument',CaseFileName);
Command=['Call getDocument().save("',CaseFileName,'")'];
invoke(MN6, 'processCommand', Command);



