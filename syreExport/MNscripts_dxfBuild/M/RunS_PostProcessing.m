%% VERSIONE 20 11 2011
% Post Processing:
% read the global quantities (torque, flux, loss...)
% Output:
%   - Ft (time Vs Flux, Current, Torque ..)
%   - save Ft Cas Sim

%% debug PostProcessing
% clear all
% close all
% [CaseFileName,FilePath] = uigetfile('*.mn','Select the ParMachine Mat-file');   %
% MN6 = actxserver('Magnet.application'); %
% set(MN6, 'Visible', 1);
% DocSV = invoke(MN6, 'openDocument',[FilePath CaseFileName]);    %
%% end debug PostProcessing

Doc         = invoke(MN6, 'getDocument');
View        = invoke(Doc, 'getCurrentView');
Solution    = invoke(Doc, 'getSolution');
calculator  = invoke(MN6, 'getStackCalculator');

% read simulation time values
NTime       = invoke(Solution,'getFieldSolutionTimeInstants',1);  % numero istanti

Command     = 'TimeInstants = getDocument().getSolution().getFieldSolutionTimeInstants(1,time_instants)';
invoke(MN6, 'processCommand', Command);

time = [];
for ij = 1:NTime
    invoke(MN6, 'processCommand', ['Call setVariant(0, time_instants(' num2str(ij-1) '))']);
    time(ij) = invoke(MN6, 'getVariant', 0);
end
%%

Nbodies = invoke(Doc,'getNumberOfBodies',1);
invoke(MN6, 'processCommand', 'REDIM ArrayOfValues(0)');

% IPM(1) o SMPM (0)
IPM = isfield(Mac,'tipo_strati');

%% 21 nov 2011
if isfield(Mac,'n_mag_simulati')
    n_mag = Mac.n_mag_simulati;
else
    if IPM
        n_mag = size(Mac.magneti,1)/3;
    else %if motor is an SPM model
        n_mag = size(Mac.magneti,1)/4;
    end
end

Flux_1sim = [];
Torque_1sim = [];
pm_loss = zeros(n_mag,NTime);

Ncoils = invoke(Doc,'getNumberOfCoils');

%% MATTEO solo per S1 ed S4 settare Ncoils %% Solo per il caso S1 S4 Matteo RICORDATELO

Ncoils=3;
Porz_stat=Mac.Qs;
Nbarrier=0;
%%
rotor_barrier_loss=zeros(NTime,Nbarrier);

%%
for ii=1:NTime
    % flux, currents, resitance
    Flux_1sim_1ph = [];
    Torque_r = [];
    
    % problem ID
    Command = 'ReDim ProblemID(1)';
    invoke(MN6, 'processCommand', Command);
    Command = 'ProblemID(0) = 1';
    invoke(MN6, 'processCommand', Command);
    Command = ['ProblemID(1) = ' num2str(time(ii) * 1.001)];
    invoke(MN6, 'processCommand', Command);
    
    for jj=1:Ncoils
        Command = ['Call getDocument.getSolution.getFluxLinkageThroughCoil(ProblemID, ' num2str(jj) ', magnitude, phase)'];
        invoke(MN6, 'processCommand', Command);
        invoke(MN6, 'processCommand', 'Call setVariant(0, magnitude)');
        x = invoke(MN6, 'getVariant', 0);
        Flux_1sim_1ph = [Flux_1sim_1ph x];
    end
    
    % torque
    for jjj=1:Nbodies
        Command = ['Call getDocument.getSolution.getTorqueAboutOriginOnBody(ProblemID, ' num2str(jjj) ', torque_x, torque_y, torque_z)'];
        invoke(MN6, 'processCommand', Command);
        invoke(MN6, 'processCommand', 'Call setVariant(0, torque_z)');
        torque = invoke(MN6, 'getVariant', 0);
        Torque_r = [Torque_r torque];
    end
    
%     % co-energy
%     Command = ['Coenergy = getDocument.getSolution.getCoenergy(ProblemID)'];
%     invoke(MN6, 'processCommand', Command);
%     invoke(MN6, 'processCommand', 'Call setVariant(0, Coenergy)');
%     CoEn(ii) = invoke(MN6, 'getVariant', 0);
%    

%% load iron and pm losses (only for 360 deg simulation)
    if (ii == 1) && (Sim.sim_angle == 360)
        
        % stator: hysteresis and classical
%         h_loss_s = zeros(1,Mac.Qs);
%         c_loss_s = zeros(1,Mac.Qs);
        
        %% prendo solo 2 cave invece di Qs perchè tutte ci mette una vita
        for jjj = 1:Mac.Qs
%             for jjj = 1:Porz_stat
                Command = ['IronLoss = getDocument.getSolution.getIronLossInComponent (ProblemID,"statore_' num2str(jjj) '" , Losses)'];
                invoke(MN6, 'processCommand', Command);
                invoke(MN6, 'processCommand', 'Call setVariant(0, Losses)');
                test0 = invoke(MN6, 'getVariant', 0);
                test0 = cell2mat(test0);
                h_loss_s(jjj) = test0(1);
                c_loss_s(jjj) = test0(2);
          end
        
        % rotor
        Command = ['IronLoss = getDocument.getSolution.getIronLossInComponent (ProblemID,"rotor" , Losses)'];
        invoke(MN6, 'processCommand', Command);
        invoke(MN6, 'processCommand', 'Call setVariant(0, Losses)');
        test0 = invoke(MN6, 'getVariant', 0);
        test0 = cell2mat(test0);
        h_loss_r = test0(1);
        c_loss_r = test0(2);
    end
    
    %% PM loss
    for jjj = 1:n_mag
        Command = ['PowerLoss=getDocument().getSolution().getOhmicLossInConductor(ProblemID,"magnet_' num2str(jjj) '")'];
        invoke(MN6, 'processCommand', Command);
        invoke(MN6, 'processCommand', 'Call setVariant(0, PowerLoss)');
        pm_loss(jjj,ii) = invoke(MN6, 'getVariant', 0);
    end
    Flux_1sim = [Flux_1sim;Flux_1sim_1ph];
    Torque_1sim = [Torque_1sim;Torque_r];
    ii;
%     keyboard
    %% Rotor conductor Losses:
    if (Nbarrier==0)
        jjj=1;
        rotor_barrier_loss(ii,jjj) = 0;
    else
%         keyboard
        jk=1;
        for jjj=1:Nbarrier
            for jj=1:Mac.ps
                Command = ['PowerLoss=getDocument().getSolution().getOhmicLossInConductor(ProblemID,"barrier#' num2str(jjj),'_',num2str(jj) '")'];
                invoke(MN6, 'processCommand', Command);
                invoke(MN6, 'processCommand', 'Call setVariant(0, PowerLoss)');
                rotor_barrier_loss(ii,jk) = invoke(MN6, 'getVariant', 0);
                jk=jk+1;
            end
        end
    end
    
    
end
% keyboard
Torque_1sim = Torque_1sim * Mac.Q /Mac.Qs;
Flux_1sim = Flux_1sim * Mac.Q /Mac.Qs / Cas.N_parallel;
% CoEnergy = CoEn * Mac.Q /Mac.Qs;
pm_loss = pm_loss * Mac.Q /Mac.Qs;
Ppm = mean(sum(pm_loss));
if (Sim.sim_angle == 360)
%     Pfes_h = mean(h_loss_s) * Mac.Q;
%     Pfes_c = mean(c_loss_s) * Mac.Q;
%     Pfer_h = h_loss_r * 2 * Mac.p;
%     Pfer_c = c_loss_r * 2 * Mac.p;
%%  Matteo modifica preliminare
Period_mot=gcd(Mac.Qs,Mac.p);
if (Porz_stat==1)
    Pfes_h = mean(h_loss_s) * Mac.Q /Mac.Qs;
    Pfes_c = mean(c_loss_s) * Mac.Q /Mac.Qs;
else
    Pfes_h = mean(h_loss_s) * Mac.Q;
    Pfes_c = mean(c_loss_s) * Mac.Q;

end
%     Pfer_h = h_loss_r * 2 * Mac.p/Period_mot;     % usati nel caso delle macchine S1 ed S2 possono essere erronee come linee di prog
%     Pfer_c = c_loss_r * 2 * Mac.p/Period_mot;
    
    Pfer_h = h_loss_r / Mac.ps * 2 * Mac.p;
    Pfer_c = c_loss_r / Mac.ps * 2 * Mac.p;
  
    PjrBar=sum(mean(rotor_barrier_loss))*2 * Mac.p/Mac.ps;

%%
end

% --> Fluxabc, Fluxdq (Npos x 1)
RunS_FluxElab
RunS_VoltageElab

Torque_0 = -mean(Torque_1sim(:,1));
% SaveData - single simulation
tmp1=Cas.Id_(k)*ones(size(time,2),1);
tmp2=Cas.Iq_(k)*ones(size(time,2),1);

% Ft.values=[time(2:end)' Flux_1sim(2:end,:) Ea Eb Ec tmp1(2:end) tmp2(2:end) Fluxd Fluxq Torque_1sim((2:end),:) CoEnergy(2:end)'];
Ft.values=[time(2:end)' Flux_1sim(2:end,:) Ea Eb Ec tmp1(2:end) tmp2(2:end) Fluxd Fluxq Torque_1sim((2:end),:)];
if (Sim.sim_angle == 360)
    save([CaseFileName(1:end-3) '.mat'],'Ft','Cas','Mac','Sim','pm_loss','Ppm','Pfes_h','Pfes_c','Pfer_h','Pfer_c','PjrBar');
else
    save([CaseFileName(1:end-3) '.mat'],'Ft','Cas','Mac','Sim','pm_loss','Ppm');
end
clear tmp1 tmp2
