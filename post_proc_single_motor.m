%% pproc_single_motor.m

% FEA evaluates one machine in single or multiple id, iq conditions

% input:
% - mot_xx.mat: mat file produced by evalParetoFront.m, describing the xx motor
% under test
% - io [p.u.], gamma [deg]

%% singt mode
% simulates single or multiple id, iq conditions
% example inputs for single point simulation:
% io = 1, gamma = 45 simulates id = io * cosd(45), iq = io * sind(45)
% example of inputs for multiple points simulation:
% io = [1 2 3], gamma = [45 45 45], simulates same phase angle and variable
% amplitude .. and so on

%% singm mode (occurs when gamma = 1000)
% simulates a regular grid of (id,iq) combinations, to construct the motor
% magnetic model (d,q flux linkages over id,iq)
% example of inputs:
% io = 1, gamma = 1000

% drawYN permits to avoid to draw the motor. E.g. if little manual modifications
% have been made on the .fem file like fillets etc ..
% drawYN = 'Y' : loads the geo from mot_xx.mat and re-draws the motor
% drawYN = 'N' : loads the geo from mot_xx.mat and the motor from
% mot_xx.fem without redrawing it

clear all; close all; clc;
addpath('d:\femm42_beta\mfiles\');
addpath mfiles
current_path=cd;

if matlabpool('Size')>0
    matlabpool close force
end
% matlabpool(4)
%%%%%%%%%%%%%%%%%%%
global eval_type

load results\lastpathname.mat;

%% input dialog box
temp = inputdlg({'Br [T]';'gamma [deg]';'current load [p.u.]';...
    'draw the motor [Y/N]'},'INPUT',1,{'0';'0';'10';'N'});

Br = eval(cell2mat(temp(1)));           % remanence of all barriers magnets
gamma_temp = eval(cell2mat(temp(2)));   % current phase angle
io = eval(cell2mat(temp(3)));  % current to be simulated
drawYN = cell2mat(temp(4));

% assign to variable geo.
% geo.Ns=Ns;
% geo.Nbob=geo.Ns/geo.p/(geo.ns/6)/size(geo.avv,1);  % conductors in slot per label

if gamma_temp == 1000
    eval_type = 'singm';
else
    eval_type = 'singt';
end

openfemm

switch drawYN
    
    case 'Y'
        
        [filemot, pathname] = uigetfile([pathname '\*m*.mat'], 'Pick a motor');
        save results\lastpathname.mat pathname;
        load([pathname filemot]);
        draw_motor_in_FEMM;
        
    case 'N'
        
        [filemot, pathname] = uigetfile([pathname '\*.fem'], 'Pick a motor');
        save results\lastpathname.mat pathname;
        copyfile([pathname filemot],[current_path '\mot0.fem']);
        filemot = strrep(filemot,'.fem','.mat');
        load([pathname filemot]);
end

% Magnet coercivity
        Hc = 1/(4e-7*pi)*Br;    % PM coercivity
        geo.Hc = Hc;

% calc current amplitude for simultation with one turns
io_femm=io*geo.Nbob;


switch eval_type
    
    case 'singt'
        
        gammat = gamma_temp;
        
        for i = 1:length(io)
            
            % run FEMM
            out = eval_motor_in_FEMM(geo,io_femm(i),gammat(i));
            
            Istr=num2str(io(i),3); Istr = strrep(Istr,'.','A');
            gammaStr=num2str(gamma_temp(i),4); gammaStr = strrep(gammaStr,'.','d');            
            if isempty(strfind(gammaStr, 'd'))
                gammaStr = [gammaStr 'd'];
            end
            FILENAME = [filemot(1:end-4) '_T_eval_',Istr,'_',gammaStr];
            [success,message,messageid] = mkdir(pathname,FILENAME);
            NewDir=[pathname,FILENAME,'\'];
            copyfile('sim_mot_temp.mat',[NewDir,FILENAME,'.mat']);
            save([NewDir,FILENAME,'.mat'],'geo','-append');
            
            % plots the waveforms and saves the figs into the new folder
            klength = 1; kturns = 1;
            plot_singt
            
        end
        
    case 'singm'
        
        n = 3;
        
        switch length(io)
            
            case 1
                idvect = linspace(0,io,n); 
                iqvect = linspace(0,io,n);
            case 2
                idvect = linspace(0,io(1),n);
                iqvect = linspace(0,io(2),n);
        
        end
         
        % run FEMM over the id,iq grid
        F_map = eval_FdFq_tables_in_FEMM(geo,idvect,iqvect,0);
        
        % builds a new folder for each id, iq simulation
        Idstr=num2str(max(abs(idvect)),3); Idstr = strrep(Idstr,'.','A');
        Iqstr=num2str(max(abs(iqvect)),3); Iqstr = strrep(Iqstr,'.','A');
        
        if isempty(strfind(Idstr, 'A'))
            Idstr = [Idstr 'A'];
        end
        if isempty(strfind(Iqstr, 'A'))
            Iqstr = [Iqstr 'A'];
        end
        FILENAME = [filemot(1:end-4) '_F_map'];
        [success,message,messageid] = mkdir(pathname,[FILENAME '_' Idstr 'x' Iqstr]);
        NewDir=[pathname,[FILENAME '_' Idstr 'x' Iqstr],'\'];
        save([NewDir,FILENAME,'.mat'],'F_map');
                    
        % interp and then plots the magnetic curves
        klength = 1; kturns = 1;
        plot_singm
        
end

closefemm
matlabpool close;
