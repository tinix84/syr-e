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
% Uses matlabpool (parfor)

%   Key INPUTs: CurrLoPP: current to be simulated
%               GammaPP: current phase angle
%               BrPP: remanence of all barriers magnets
%               NumOfRotPosPP: # simulated positions
%               AngularSpanPP: angular span of simulation
%               NumGrid: number of points in [0 Imax] for the single machine post-processing
%=========================================================================

% singt mode
% simulates single or multiple (id,iq) conditions
% example inputs:
% single condition: CurrLoPP = 1, GammaPP = 45
% multiple points:  CurrLoPP = [1 1.5 2], gamma = [45 45 45]

% singm mode (occurs when gamma = 1000)
% regular grid of (id,iq) combinations, builds the magnetic model
% (d,q flux linkages over id,iq)
% example of inputs:
% CurrLoPP = 1, GammaPP = 1000

if (nargin) == 0
    % input dialog box, script run out of GUI
    [filemot, pathname, fltidx]=uigetfile(' *.fem', 'Pick a motor');
    load(strrep(filemot,'.fem','.mat'));
    if isoctave()
        temp = inputdlg({'current load [p.u.]';'gamma [deg]';'Br [T]'; ...
            'number of rotor positions';'angular span (elt. deg.)'; ...
            'points in [0 Imax]';'Evaluate loss? [1 - Yes, 0 - No]'},'INPUT',...
            1,{dataSet.CurrLoPP;dataSet.GammaPP;dataSet.BrPP;dataSet.NumOfRotPosPP;dataSet.AngularSpanPP;dataSet.NumGrid;dataSet.LossEvaluationCheck});
    else
        temp = inputdlg({'current load [p.u.]';'gamma [deg]';'Br [T]'; ...
            'number of rotor positions';'angular span (elt. deg.)'; ...
            'points in [0 Imax]';'Evaluate loss? [1 - Yes, 0 - No]'},'INPUT', ...
            1,{num2str(dataSet.CurrLoPP);num2str(dataSet.GammaPP);num2str(dataSet.BrPP);num2str(dataSet.NumOfRotPosPP);num2str(dataSet.AngularSpanPP);num2str(dataSet.NumGrid);num2str(dataSet.LossEvaluationCheck)});
    end
    dataIn.LossEvaluationCheck = eval(cell2mat(temp(7)));
    
    CurrLoPP = eval(cell2mat(temp(1)));         % current to be simulated
    GammaPP = eval(cell2mat(temp(2)));          % current phase angle
    BrPP = eval(cell2mat(temp(3)));             % remanence of all barriers magnets
    NumOfRotPosPP = eval(cell2mat(temp(4)));    % # simulated positions
    AngularSpanPP = eval(cell2mat(temp(5)));    % angular span of simulation
    NumGrid = eval(cell2mat(temp(6)));          % number of points in [0 Imax] for the single machine post-processing
    
    % Iron Loss Input
    if dataIn.LossEvaluationCheck == 1
        temp = inputdlg({'Evaluated speed [rpm]'},'INPUT',1,{'0'});
        per.EvalSpeed = eval(cell2mat(temp(1)));
    end
    
else
    
    pathname=dataIn.currentpathname;
    filemot = strrep(dataIn.currentfilename,'.mat','.fem');
    load([dataIn.currentpathname dataIn.currentfilename]);
    
    CurrLoPP = dataIn.CurrLoPP;
    GammaPP  = dataIn.GammaPP;
    BrPP = dataIn.BrPP;
    NumOfRotPosPP = dataIn.NumOfRotPosPP;
    AngularSpanPP = dataIn.AngularSpanPP;
    NumGrid = dataIn.NumGrid;
    
    % Iron Loss Input
    if dataIn.LossEvaluationCheck == 1
        per.EvalSpeed = dataIn.EvalSpeed;
    end
    
end

clc;

% overload_temp =  CurrLoPP;   % current to be simulated
gamma_temp = GammaPP;        % current phase angle
% Br = BrPP;                   % remanence of all barriers magnets

% if gamma_temp == 1000
if GammaPP == 1000
    eval_type = 'singm';    % map over id, iq
else
    eval_type = 'singt';    % single or multiple id, iq combination
end

per.overload=CurrLoPP;
per.BrPP=BrPP;

geo.nsim_singt = NumOfRotPosPP;       % # simulated positions
geo.delta_sim_singt = AngularSpanPP;  % angular span of simulation

iAmp = dataIn.SimulatedCurrent;

switch eval_type
    
    case 'singt'
        
        % single point or array of points simulation
        performance = cell(1,length(CurrLoPP));
        output = cell(1,length(CurrLoPP));
        geometry = cell(1,length(CurrLoPP));
        tempDirName = cell(1,length(CurrLoPP));
        for i = 1:length(CurrLoPP)
            performance{i} = per;
            performance{i}.overload = CurrLoPP(i);
%             performance{i}.gamma=gamma_temp(i);
            performance{i}.gamma=GammaPP(i);
        end
        
        geo.RemoveTMPfile = 'OFF';
        for i = 1:length(CurrLoPP)
            %             if dataIn.LossEvaluationCheck == 0
            [~,geometry{i},mat,output{i},tempDirName{i}] = FEMMfitness([],geo,performance{i},mat,eval_type,[pathname,filemot]);
            %             else
            %                 [~,geometry{i},mat,output{i},tempDirName{i}] = FEMMfitness_IronLoss([],geo,performance{i},mat,eval_type,[pathname,filemot]);
            %             end
        end
        
        % save output into individual folders
        %         iAmp = CurrLoPP*calc_io(geo,per);
        for i = 1:length(CurrLoPP)
            
            geo = geometry{i};
            out = output{i};
            per = performance{i};
            dirName = tempDirName{i};
            
            iStr=num2str(iAmp(i),3); iStr = strrep(iStr,'.','A');
            gammaStr=num2str(GammaPP(i),4); gammaStr = strrep(gammaStr,'.','d');
            if isempty(strfind(gammaStr, 'd'))
                gammaStr = [gammaStr 'd'];
            end
            
            FILENAME = [filemot(1:end-4) '_T_eval_',iStr,'_',gammaStr];
            [success,message,messageid] = mkdir(pathname,FILENAME);
            newDir=[pathname,FILENAME,'\'];
            
            if isoctave()            %OCT
                file_name1= strcat(newDir,FILENAME,'.mat');
                save('-mat7-binary', file_name1,'geo','per','out');
                dirIn=strcat(dirName, '\mot_temp.fem');
                dirDest=strcat(newDir, FILENAME, '.fem');
                movefile(dirIn, dirDest);
                clear file_name1 dirIn dirDest
            else
                save([newDir,FILENAME,'.mat'],'geo','per','out');
                movefile([dirName 'mot_temp.fem'],[newDir FILENAME '.fem']);
            end
            
            % lot and save figs
            klength = 1; kturns = 1; delta_sim_singt = geo.delta_sim_singt;
            plot_singt(out,klength,kturns,delta_sim_singt,newDir,filemot);
            
        end
        
        % extra figs, if input current is array
        if length(CurrLoPP)>1
            
            id = zeros(1,length(CurrLoPP));
            iq = zeros(1,length(CurrLoPP));
            T = zeros(1,length(CurrLoPP));
            dTpu = zeros(1,length(CurrLoPP));
            dTpp = zeros(1,length(CurrLoPP));
            fd = zeros(1,length(CurrLoPP));
            fq = zeros(1,length(CurrLoPP));
            
            for i = 1:length(CurrLoPP)
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
            
            x = 1:length(CurrLoPP);
            figure();
            if ~isoctave()
                figSetting();
            end
            subplot(2,1,1)
            plot(x,T,'-x',x,T+0.5*dTpp,'r',x,T-0.5*dTpp,'r'), grid on, ylabel('$T$ [Nm]')
            subplot(2,1,2)
            plot(x,dTpp,'-x'), grid on, ylabel('$\Delta T_{pp}$ [Nm]')
            xlabel('simulation \#')
            h=gcf();
            if isoctave() %OCT
                fig_name=strcat(dirPower, filemot(1:end-4), '_torque_sens');
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[dirPower,filemot(1:end-4),'_torque_sens.fig'])
            end
            
            figure()
            if ~isoctave()
                figSetting();
            end
            subplot(2,1,1)
            plot(x,fd,'-x',x,fq,'-x'), grid on, ylabel('[Vs]'), legend('$\lambda_d$','$\lambda_q$'),
            subplot(2,1,2)
            plot(x,abs(sin(atan(iq./id)-atan(fq./fd))),'-x'), grid on, ylabel('$cos \varphi$')
            xlabel('simulation \#'),
            h=gcf();
            if isoctave() %OCT
                fig_name=strcat(dirPower, filemot(1:end-4), '_fdq_IPF_sens');
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[dirPower,filemot(1:end-4),'_fdq_IPF_sens.fig'])
            end
            
        end
        
    case 'singm'
        
        % flux map over a rectangular grid of (id,iq) combinations
        n_grid = NumGrid;     % number of points in [0 Imax] for the single machine post-processing
        %         iAmp = overload_temp*calc_io(geo,per);
        
        switch length(iAmp)
            case 1  % square domain
                if (strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'Vtype'))
                    idvect = linspace(-iAmp,0,n_grid);
                    iqvect = linspace(0,iAmp,n_grid);
                else
                    idvect = linspace(0,iAmp,n_grid);
                    iqvect = linspace(0,iAmp,n_grid);
                end
            case 2  % rectangular domain
                if (strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'Vtype'))
                    idvect = linspace(-iAmp(1),0,n_grid);
                    iqvect = linspace(0,iAmp(2),n_grid);
                else
                    idvect = linspace(0,iAmp(1),n_grid);
                    iqvect = linspace(0,iAmp(2),n_grid);
                end
            case 4  % CurrLoPP = [IdMin IdMax IqMin IqMax]
                idvect=linspace(iAmp(1),iAmp(2),n_grid);
                iqvect=linspace(iAmp(3),iAmp(4),n_grid);
        end
        
        [F_map,OUT] = eval_FdFq_tables_in_FEMM(geo,per,idvect,iqvect,eval_type,[pathname filemot],mat,dataIn);
        
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
        clear pathname; %AS
        if isoctave()            %OCT
            file_name1= strcat(NewDir,FILENAME,'.mat');
            save('-v7', file_name1,'F_map','OUT');
            clear file_name1
        else
            save([NewDir,FILENAME,'.mat'],'F_map','OUT');
        end
        
        % interp and then plots the magnetic curves
        n_interp = 256;    % number of points in [0 Imax] for data interpolation
        klength = 1; kturns = 1; n2=n_interp;
        filemot=filemot(1:end-4);
        plot_singm;
        
end


