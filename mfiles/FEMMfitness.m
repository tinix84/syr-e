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

function [cost,geo,mat,out,pathname] = FEMMfitness(RQ,geo,per,mat,eval_type,filename)

[thisfilepath,pathname]=createTempDir();

if ~isempty(RQ)
    
    % MODE optimization (RQ geometry)
    RQ=RQ';
    
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
    if ~strcmp(eval_type,'singt')
        RQ
    end
    % debug .. end
    [geo,gamma,mat] = interpretRQ(RQ,geo,mat);
    per.gamma=gamma;
    
    openfemm
    [geo,mat] = draw_motor_in_FEMM(geo,eval_type,mat);
    mi_saveas([pathname 'mot0.fem']);
else
    % post proc or FEMM simulation (existing geometry)
    openfemm
    opendocument(filename);
    mi_saveas([pathname 'mot0.fem']); % saves the file with name ’filename’
end

% evaluates the candidate machine (T, DT, fd, fq)
fem = dimMesh(geo,eval_type);
filename='mot0.fem';
mat.LayerMag.Br = per.BrPP;
mat.LayerMag.Hc = per.BrPP/(4e-7*pi*mat.LayerMag.mu);

[SOL] = simulate_xdeg(geo,per,mat,eval_type,pathname,filename,0);

out.id = mean(SOL.id);
out.iq = mean(SOL.iq);
out.fd = mean(SOL.fd);
out.fq = mean(SOL.fq);
out.T= abs(mean(SOL.T));
out.dT = std(SOL.T);
out.dTpu = std(SOL.T)/out.T;
out.dTpp = max(SOL.T)-min(SOL.T);
out.IPF = sin(atan(out.iq./out.id)-atan(out.fq./out.fd));
out.MassPM=mean(SOL.VolPM)*mat.LayerMag.kgm3; %[kg] massa magneti rotore - rev.Gallo 20/03/2018
out.SOL = SOL;

if isfield(SOL,'F')
    out.F=mean(SOL.F);
end

if isfield(SOL,'Pfes_h')
    out.Pfes_h = SOL.Pfes_h;
    out.Pfes_c = SOL.Pfes_c;
    out.Pfer_h = SOL.Pfer_h;
    out.Pfer_c = SOL.Pfer_c;
    out.Pfe_total = out.Pfes_h + out.Pfes_c + out.Pfer_h + out.Pfer_c;
end

if ~isempty(RQ)
    
    % Cost functions
    cost = zeros(1,length(geo.OBJnames));
    temp1 = 1;
    % Torque
    if strcmp(geo.OBJnames{temp1},'Torque')
        cost(temp1) = -out.T;
        temp1 = temp1+1;
    end
    % Torque Ripple
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'TorRip')
        %         cost(temp1) = out.dTpu*100;
        cost(temp1) = out.dTpp;
        temp1 = temp1+1;
    end
    % Copper Mass
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'MassCu')
        cost(temp1)=calcMassCu(geo);
        temp1=temp1+1;
    end
    % PM Mass - rev.Gallo 20/03/2018
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'MassPM')
        cost(temp1) = out.MassPM;
    end
    
    % penalize weak solutions
    for j = 1:length(cost)
        if cost(j)>per.objs(j,1)
            if per.objs(j,1)>0
                cost(j)=cost(j)*10;  % minimization problem
            else
                cost(j)=cost(j)*0.1; % maximization problem
            end
        end
    end
    
    %     dataSetPath = strcat(thisfilepath,'\dataSet.mat');    %OCT
    load('dataSet.mat');
    geo.RQ = RQ;
    
    [dataSet] = SaveInformation(geo,mat,dataSet);
    %delete('dataSet.mat');
    if isoctave()            %OCT
        save('-mat7-binary', 'mot0.mat','geo','cost','per','dataSet','mat');
    else
        save([pathname 'mot0.mat'],'geo','cost','per','dataSet','mat');
    end
else
    cost = [];
    %save([pathname 'geo'],'geo','out','mat','-append');   % save geo and out
    save([pathname 'geo.mat'],'geo','out','mat');   % save geo and out
end

% cd(currentDir);
