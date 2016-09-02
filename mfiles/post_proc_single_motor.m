% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function post_proc_single_motor(dataIn,varargin)

% post_proc_single_motor simulates an existing machine
% Open matlabpool manually prior to execution

%% DATA ===================================================================
%%   INPUT: CurrLoPP: current to be simulated
%%          GammaPP: current phase angle
%%          BrPP: remanence of all barriers magnets
%%          NumOfRotPosPP: # simulated positions
%%          AngularSpanPP: angular span of simulation
%%          NumGrid: number of points in [0 Imax] for the single machine post-processing
%%=========================================================================

%% singt mode
% simulates single or multiple (id,iq) conditions
% example inputs:
% single condition: CurrLoPP = 1, GammaPP = 45
% multiple points:  CurrLoPP = [1 1.5 2], gamma = [45 45 45]

%% singm mode (occurs when gamma = 1000)
% regular grid of (id,iq) combinations, builds the magnetic model
% (d,q flux linkages over id,iq)
% example of inputs:
% CurrLoPP = 1, GammaPP = 1000

%% input dialog box
if (nargin) == 0
    %% input dialog box
    temp = inputdlg({'current load [p.u.]';'gamma [deg]';'Br [T]';'number of rotor positions';'angular span (elt. deg.)';'points in [0 Imax]'},'INPUT',1,{'[1 1 1 1 1 1 1 1]';'[40 45 50 55 60 65 70 75]';'0';'6';'60';'5'});
    
    CurrLoPP = eval(cell2mat(temp(1)));         % current to be simulated
    GammaPP = eval(cell2mat(temp(2)));          % current phase angle
    BrPP = eval(cell2mat(temp(3)));             % remanence of all barriers magnets
    NumOfRotPosPP = eval(cell2mat(temp(4)));    % # simulated positions
    AngularSpanPP = eval(cell2mat(temp(5)));    % angular span of simulation
    NumGrid = eval(cell2mat(temp(6)));          % number of points in [0 Imax] for the single machine post-processing
else
    CurrLoPP = dataIn.CurrLoPP;
    GammaPP  = dataIn.GammaPP;
    BrPP = dataIn.BrPP;
    NumOfRotPosPP = dataIn.NumOfRotPosPP;
    AngularSpanPP = dataIn.AngularSpanPP;
    NumGrid = dataIn.NumGrid;
end

clc;
syreRoot = fileparts(which('GUI_Syre.m'));
current_path = syreRoot;

overload_temp =  CurrLoPP;   % current to be simulated
gamma_temp = GammaPP;        % current phase angle
Br = BrPP;                   % remanence of all barriers magnets

if gamma_temp == 1000
    eval_type = 'singm';    % map over id, iq
else
    eval_type = 'singt';    % single or multiple id, iq combination
end

% choose the .fem file and load the associated .mat file
% [filemot, pathname] = uigetfile([syreRoot '\*.fem'], 'Pick a motor');
pathname = dataIn.currentpathname;
filemot = strrep(dataIn.currentfilename,'.mat','.fem');
load([pathname strrep(filemot,'.fem','.mat')]);

geo.nsim_singt = NumOfRotPosPP;       % # simulated positions
geo.delta_sim_singt = AngularSpanPP;  % angular span of simulation

%% Iron Loss Input
if dataIn.LossEvaluationCheck == 1
%     geo.loss.kh = dataIn.HysteresisLossFactor;
%     geo.loss.alpha = dataIn.HysteresisFrequencyFactor;
%     geo.loss.beta = dataIn.HysteresisFluxDenFactor;
%     geo.loss.ke = dataIn.EddyCurLossFactor;
%     geo.rhoFE = dataIn.IronMassDen;
    per.EvalSpeed = dataIn.EvalSpeed;
end

% Magnet coercivity
Hc = 1/(4e-7*pi)*Br;    % PM coercivity
mat.LayerMag.Hc = Hc;

switch eval_type
    
    case 'singt'
        
        % single or multiple simulation
        performance = cell(1,length(overload_temp));
        output = cell(1,length(overload_temp));
        geometry = cell(1,length(overload_temp));
        tempDirName = cell(1,length(overload_temp));
        for i = 1:length(overload_temp)
            performance{i} = per;
            performance{i}.overload = overload_temp(i);
            performance{i}.gamma=gamma_temp(i);
        end
        
        geo.RemoveTMPfile = 'OFF';
        for i = 1:length(overload_temp)
            if dataIn.LossEvaluationCheck == 0
                [~,geometry{i},mat,output{i},tempDirName{i}] = FEMMfitness([],geo,performance{i},mat,eval_type,[pathname,filemot]);
            else
                [~,geometry{i},mat,output{i},tempDirName{i}] = FEMMfitness_IronLoss([],geo,performance{i},mat,eval_type,[pathname,filemot]);
            end
        end
        
        % save output into individual folders
        iAmp = overload_temp*calc_io(geo,per);
        for i = 1:length(overload_temp)
            
            geo = geometry{i};
            out = output{i};
            per = performance{i};
            dirName = tempDirName{i};
            
            iStr=num2str(iAmp(i),3); iStr = strrep(iStr,'.','A');
            gammaStr=num2str(gamma_temp(i),4); gammaStr = strrep(gammaStr,'.','d');
            if isempty(strfind(gammaStr, 'd'))
                gammaStr = [gammaStr 'd'];
            end
            
            FILENAME = [filemot(1:end-4) '_T_eval_',iStr,'_',gammaStr];
            [success,message,messageid] = mkdir(pathname,FILENAME);
            newDir=[pathname,FILENAME,'\'];
            
            save([newDir,FILENAME,'.mat'],'geo','per','out');
            movefile([syreRoot '\tmp\' dirName '\mot_temp.fem'],[newDir FILENAME '.fem']);
            %             plot and save figs
            klength = 1; kturns = 1; delta_sim_singt = geo.delta_sim_singt;
            plot_singt(out,klength,kturns,delta_sim_singt,newDir,filemot);
        end
        
        % extra figs, if input current is array
        if length(overload_temp)>1
            
            id = zeros(1,length(overload_temp));
            iq = zeros(1,length(overload_temp));
            T = zeros(1,length(overload_temp));
            dTpu = zeros(1,length(overload_temp));
            dTpp = zeros(1,length(overload_temp));
            fd = zeros(1,length(overload_temp));
            fq = zeros(1,length(overload_temp));
            
            for i = 1:length(overload_temp)
                id(i) = output{i}.id;
                iq(i) = output{i}.iq;
                T(i) = output{i}.T;
                dTpu(i) = output{i}.dTpu;
                dTpp(i) = output{i}.dTpp;
                fd(i) = output{i}.fd;
                fq(i) = output{i}.fq;
            end
            dirPower=[pathname,filemot(1:end-4),'singT\'];
            mkdir(dirPower);
            
            x = 1:length(overload_temp);
            figure(10), subplot(2,1,1)
            plot(x,T,'-x',x,T+0.5*Tpp,'r',x,T-0.5*dTpp,'r'), grid on, ylabel('torque [Nm]')
            subplot(2,1,2)
            plot(x,dTpp,'-x'), grid on, ylabel('torque ripple pk-pk [Nm]')
            xlabel('simulation #'),
            saveas(gcf,[dirPower,filemot(1:end-4),'_torque_sens.fig'])
            
            figure(11), subplot(2,1,1)
            plot(x,fd,'-x',x,fq,'-x'), grid on, ylabel('[Vs]'), legend('\lambda_d','\lambda_q'),
            subplot(2,1,2)
            plot(x,sin(atan(iq./id)-atan(fq./fd)),'-x'), grid on, ylabel('IPF')
            xlabel('simulation #'),
            saveas(gcf,[dirPower,filemot(1:end-4),'_fdq_IPF_sens.fig'])
            
        end
        
    case 'singm'
        
        % flux map over a grid of id,iq combinations
        n_grid = NumGrid;     % number of points in [0 Imax] for the single machine post-processing
        iAmp = overload_temp*calc_io(geo,per);
        switch length(iAmp)
            case 1  % square domain
                if strcmp(geo.RotType,'SPM')
                    idvect = linspace(-iAmp,0,n_grid);
                    iqvect = linspace(0,iAmp,n_grid);
                else
                    idvect = linspace(0,iAmp,n_grid);
                    iqvect = linspace(0,iAmp,n_grid);
                end
            case 2  % rectangular domain
                if strcmp(geo.RotType,'SPM')
                    idvect = linspace(-iAmp(1),0,n_grid);
                    iqvect = linspace(0,iAmp(2),n_grid);
                else
                    idvect = linspace(0,iAmp,n_grid);
                    iqvect = linspace(0,iAmp,n_grid);
                end
        end
        
        F_map = eval_FdFq_tables_in_FEMM(geo,per,idvect,iqvect,eval_type,[pathname filemot]);
        
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
        n_interp = 256;    % number of points in [0 Imax] for data interpolation
        klength = 1; kturns = 1; n2=n_interp;
        filemot=filemot(1:end-4);
        plot_singm;
        
end


