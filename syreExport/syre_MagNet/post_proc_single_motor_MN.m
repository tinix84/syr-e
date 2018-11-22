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

function post_proc_single_motor_MN(dataIn,varargin)

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
    [filemot, pathname, fltidx]=uigetfile(' *.fem', 'Pick a motor');
    load(strrep(filemot,'.fem','.mat'));
    if isoctave()
        temp = inputdlg({'current load [p.u.]';'gamma [deg]';'Br [T]';'number of rotor positions';'angular span (elt. deg.)';'points in [0 Imax]';'Evaluate loss? [1 - Yes, 0 - No]'},'INPUT',1,{dataSet.CurrLoPP;dataSet.GammaPP;dataSet.BrPP;dataSet.NumOfRotPosPP;dataSet.AngularSpanPP;dataSet.NumGrid;dataSet.LossEvaluationCheck});
    else
        temp = inputdlg({'current load [p.u.]';'gamma [deg]';'Br [T]';'number of rotor positions';'angular span (elt. deg.)';'points in [0 Imax]';'Evaluate loss? [1 - Yes, 0 - No]'},'INPUT',1,{num2str(dataSet.CurrLoPP);num2str(dataSet.GammaPP);num2str(dataSet.BrPP);num2str(dataSet.NumOfRotPosPP);num2str(dataSet.AngularSpanPP);num2str(dataSet.NumGrid);num2str(dataSet.LossEvaluationCheck)});
    end
    dataIn.LossEvaluationCheck = eval(cell2mat(temp(7)));
    
    CurrLoPP = eval(cell2mat(temp(1)));         % current to be simulated
    GammaPP = eval(cell2mat(temp(2)));          % current phase angle
    BrPP = eval(cell2mat(temp(3)));             % remanence of all barriers magnets
    NumOfRotPosPP = eval(cell2mat(temp(4)));    % # simulated positions
    AngularSpanPP = eval(cell2mat(temp(5)));    % angular span of simulation
    NumGrid = eval(cell2mat(temp(6)));          % number of points in [0 Imax] for the single machine post-processing
    
    %[dataSet,~,~] = back_compatibility(dataSet,geo,per);
    %name_file1 = strrep(filemot,'.fem','.mat');
    %if isoctave()  %OCT
    %    name_file = strcat(pathname, name_file1);
    %    save ('-mat7-binary', name_file,'geo','per','dataSet','mat');
    %else
    %    save([pathname name_file1],'geo','per','dataSet','mat');
    %end
    %clear name_file1
else
    
    pathname=dataIn.currentpathname;
    %  filemot = strrep(dataIn.currentfilename,'.mat','.fem');
    filemot= dataIn.currentfilename;
    load([dataIn.currentpathname dataIn.currentfilename]);
        
    CurrLoPP = dataIn.CurrLoPP;
    GammaPP  = dataIn.GammaPP;
    BrPP = dataIn.BrPP;
    NumOfRotPosPP = dataIn.NumOfRotPosPP;
    AngularSpanPP = dataIn.AngularSpanPP;
    NumGrid = dataIn.NumGrid;
    per.EvalSpeed = dataIn.EvalSpeed;
    
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

per.overload=CurrLoPP;
per.BrPP=BrPP;

geo.nsim_singt = NumOfRotPosPP;       % # simulated positions
geo.delta_sim_singt = AngularSpanPP;  % angular span of simulation

% Magnet coercivity (QUESTO NON SERVE IN MAGNET - CONTANO SOLO LA TEMPERATURA DEL MAGNETE E LA DEFINIZIONE DEL MATERIALE)
% Hc = 1/(4e-7*pi)*Br;    % PM coercivity
% mat.LayerMag.Hc = Hc;

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
            [geometry{i},mat,output{i},tempDirName{i}] = MNfitness([],geo,performance{i},mat,eval_type,pathname,filemot);
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
            
            FILENAME = [filemot(1:end-4) '_T_eval_',iStr,'_',gammaStr '_MN'];
            %             [success,message,messageid] = mkdir(pathname,FILENAME);
            newDir=[pathname,FILENAME,'\'];
            
            %             if isoctave()            %OCT
            %                 file_name1= strcat(newDir,FILENAME,'.mat');
            %                 save('-mat7-binary', file_name1,'geo','per','out');
            %                 dirIn=strcat(dirName, '\mot_temp.fem');
            %                 dirDest=strcat(newDir, FILENAME, '.fem');
            %                 movefile(dirIn, dirDest);
            %                 clear file_name1 dirIn dirDest
            %             else
            save([newDir,FILENAME,'.mat'],'geo','per','out');
            %                 %movefile([syreRoot '\tmp\' dirName '\mot_temp.fem'],[newDir FILENAME '.fem']);
            %                 movefile([dirName 'mot_temp.fem'],[newDir FILENAME '.fem']);
            %             end
            
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
            plot(x,T,'-x',x,T+0.5*dTpp,'r',x,T-0.5*dTpp,'r'), grid on, ylabel('torque [Nm]')
            subplot(2,1,2)
            plot(x,dTpp,'-x'), grid on, ylabel('torque ripple pk-pk [Nm]')
            xlabel('simulation #'),
            h=gcf();
            if isoctave() %OCT
                fig_name=strcat(dirPower, filemot(1:end-4), '_torque_sens');
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[dirPower,filemot(1:end-4),'_torque_sens.fig'])
            end
            
            figure(11), subplot(2,1,1)
            plot(x,fd,'-x',x,fq,'-x'), grid on, ylabel('[Vs]'), legend('\lambda_d','\lambda_q'),
            subplot(2,1,2)
            plot(x,sin(atan(iq./id)-atan(fq./fd)),'-x'), grid on, ylabel('IPF')
            xlabel('simulation #'),
            h=gcf();
            if isoctave() %OCT
                fig_name=strcat(dirPower, filemot(1:end-4), '_fdq_IPF_sens');
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[dirPower,filemot(1:end-4),'_fdq_IPF_sens.fig'])
            end
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
                    idvect = linspace(0,iAmp(1),n_grid);
                    iqvect = linspace(0,iAmp(2),n_grid);
                end
            case 4 % CurrLoPP = [IdMin IdMax IqMin IqMax]
                idvect=linspace(iAmp(1),iAmp(2),n_grid);
                iqvect=linspace(iAmp(3),iAmp(4),n_grid);
        end
        
        [F_map,OUT] = eval_FdFq_tables_in_MN(geo,per,idvect,iqvect,eval_type,pathname, filemot,mat,dataIn);
        
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
        %       [success,message,messageid] = mkdir(pathname,[FILENAME '_' Idstr 'x' Iqstr]);
        NewDir=[pathname,[FILENAME '_' Idstr 'x' Iqstr  '_MN'],'\'];
        clear pathname; %AS
        if isoctave()            %OCT
            file_name1= strcat(NewDir,FILENAME,'.mat');
            save('-v7', file_name1,'F_map');
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


