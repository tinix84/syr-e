% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function [cost,geo,mat,out,pathname] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot)

currentDir=pwd();
[thisfilepath,pathname]=createTempDir();

if ~isempty(RQ)
    
    % MODE optimization (RQ geometry)
    RQ=RQ';
    if  strcmp(eval_type,'MO_GA')
        RQ=RQ';
    end
    
    geo.pathname=pwd();
    
    options.iteration=0;
    options.currentgen=1;
    options.PopulationSize=1;
    
    if strcmp(eval_type,'MO_OA')
        options.iteration=options.iteration+1;
        iteration=options.iteration;
        populationSize=options.PopulationSize;
        generation=floor(iteration/populationSize)+1;
        options.currentgen=generation;
    end
    
    % debug .. when syre crashes it is useful to have visibility of last RQ
    RQ
    % debug .. end
    [geo,gamma,mat] = interpretRQ(RQ,geo,mat);
    per.gamma=gamma;
        
    openfemm
    [geo,mat] = draw_motor_in_FEMM(geo,eval_type,mat);
    
else
    % post proc or FEMM simulation (existing geometry)
    openfemm
    opendocument(filemot);
    mi_saveas('mot0.fem'); % saves the file with name ’filename’
end

% evaluates the candidate machine (T, DT, fd, fq)
fem = dimMesh(geo,eval_type);           
%[geo, ~, ~] = STATmatr(geo,fem);  
%iAmp = per.overload*calc_io(geo,per);
%iAmpCoil = iAmp*geo.Nbob*geo.n3phase; %AS

if isoctave()            %OCT save geo before sim
    save ('-mat7-binary', 'geo.mat');
else
    save 'geo.mat';
end
%[SOL] = simulate_xdeg(geo,iAmpCoil,per.BrPP,gamma,eval_type);
filename='mot0.fem';
[SOL] = simulate_xdeg(geo,per,eval_type,pathname,filename);

out.id = mean(SOL.id);
out.iq = mean(SOL.iq);
out.fd = mean(SOL.fd);
out.fq = mean(SOL.fq);
out.T= abs(mean(SOL.T));
out.dT = std(SOL.T);
out.dTpu = std(SOL.T)/out.T;
out.dTpp = max(SOL.T)-min(SOL.T);
out.IPF = sin(atan(out.iq./out.id)-atan(out.fq./out.fd));
% Variabile calcolo Massa magneti rotore - rev.Gallo 20/03/2018
out.MassPM=mean(SOL.VolPM)*mat.LayerMag.kgm3; %[kg]
out.SOL = SOL;


if isfield(SOL,'F')
    out.F=mean(SOL.F);
end

if ~isempty(RQ)
    % Cost functions
    cost = zeros(1,length(geo.OBJnames));
    temp1 = 1; %temp2 = 1;
    if strcmp(geo.OBJnames{temp1},'Torque')
        cost(temp1) = -out.T;
        temp1 = temp1+1;
    end
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'TorRip')
%         cost(temp1) = out.dTpu*100;
        cost(temp1) = out.dTpp;
        temp1 = temp1+1;
    end
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'MassCu')
        cost(temp1)=calcMassCu(geo);
        temp1=temp1+1;
    end
% MassPM inserito come funzione obiettivo da minimizzare - rev.Gallo 20/03/2018   
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'MassPM')
        cost(temp1) = out.MassPM;
    end

    % penalize weak solutions
    for j = 1:length(cost)
        %     if (cost(2)>per.max_exp_ripple || cost(1)>-per.min_exp_torque)
        if cost(j)>per.objs(j,1)
            if per.objs(j,1)>0
                cost(j)=cost(j)*10;  % minimization problem
            else
                cost(j)=cost(j)*0.1; % maximization problem
            end
        end
    end
    % end
    
    % if ~isempty(RQ)
    dataSetPath = strcat(thisfilepath,'\dataSet.mat');    %OCT
    copyfile(dataSetPath, '.');
    load('dataSet.mat');
    geo.RQ = RQ;

    [dataSet] = SaveInformation(geo,mat,dataSet);
    delete('dataSet.mat');
    if isoctave()            %OCT
        save('-mat7-binary', 'mot0.mat','geo','cost','per','dataSet','mat');
    else
        save('mot0','geo','cost','per','dataSet','mat');
    end
else
    cost = [];
    save('geo','geo','out','mat','-append');   % save geo and out
end

cd(currentDir);
