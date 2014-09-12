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

function [cost,geo,out,dirName] = FEMMfitnessX(RQ,geo,per,eval_type,filemot)

currentDir=pwd;
[thisfilepath,dirName]=createTempDir();

% if in optimization, RQ defines the candidate machine
if ~isempty(RQ)
    
    RQ=RQ';
    if  strcmp(eval_type,'MO_GA')
        RQ=RQ';
    end
    
    geo.pathname=cd;
    
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
    
    first_index = 2;
    last_index = first_index + geo.nlay - 1;
    
    dalpha_pu = RQ(first_index:last_index);
    
    % the sum of the pu angles is rescaled to one
    if sum(dalpha_pu) > 1
        dalpha_pu = dalpha_pu/sum(dalpha_pu);
    end
    
    % dalpha(2) to dalpha(nlay) in degrees
    dalpha_temp = dalpha_pu * (90/geo.p - RQ(1));
    % all dalpha in mec degrees
    geo.dalpha = [RQ(1) dalpha_temp(1:end-1)];
    
    % SPESSORE DELLE BARRIERE: 'hc_pu'
    first_index = last_index + 1;
    last_index = first_index + geo.nlay - 1;
    size(RQ);
    getComputerName();
    geo.hc_pu = RQ(first_index:last_index);
    
    if (strcmp(geo.RotType,'Fluid')||strcmp(geo.RotType,'Seg'))
        first_index = last_index + 1;
        last_index = first_index + geo.nlay - 1;
        geo.Dfe=RQ(first_index:last_index);
        
        %     first_index=last_index+1;
        %     geo.x=RQ(first_index);
        %     first_index=first_index+1;
        %     geo.kt=RQ(first_index);
        %     first_index=first_index+1;
        %     geo.b=RQ(first_index);
        %     first_index=first_index+1;
        %     geo.acs=RQ(first_index);
        %
        %     geo.xr=geo.x*geo.r;
        %     geo.lt=geo.r*(1-geo.x-geo.g/geo.r-geo.b*geo.x/p);
    end
    % current phase angle
    gamma = RQ(end);
    
    if exist([thisfilepath filesep 'empty_case.fem'],'file')>1
        FemmProblem = loadfemmfile([thisfilepath filesep 'empty_case.fem']);
    else
        FemmProblem = loadfemmfile(['c:' filesep 'empty_case.fem']);%TODO: fix this with something more robust and crossplatform
    end
    FemmProblem.ProbInfo.Depth = geo.l;
    FemmProblem.Segments= [];
    FemmProblem.ArcSegments= [];
    FemmProblem.Nodes= [];
    FemmProblem.BoundaryProps= [];
    FemmProblem.Circuits= [];
    FemmProblem.BlockLabels= [];
    FemmProblem.PointProps= [];
    % FemmProblem = initMaterials(FemmProblem, geo);
    
    [geo,FemmProblem] = draw_motor_in_XFEMM(geo,eval_type,FemmProblem);

else
    
    FemmProblem = loadfemmfile(filemot);
    copyfile([filemot],[dirName 'mot0X.fem']);

end

%%%% uncomment for graphical output
% plotfemmproblem(FemmProblem),drawnow;
%%%% uncomment for graphical output

%% evaluates the candidate machine (T, DT, fd, fq)
iAmp = per.overload*calc_io(geo,per);
if isempty(RQ)
    gamma=per.gamma;
end
iAmpCoil=iAmp*geo.Nbob;

[SOL,FemmProblem] = simulate_xdegX(FemmProblem,geo,iAmpCoil,gamma,eval_type);

out.id = mean(SOL(:,2));
out.iq = mean(SOL(:,3));
out.fd = mean(SOL(:,4));
out.fq = mean(SOL(:,5));
out.T= abs(mean(SOL(:,6)));
out.dTpu = std(SOL(:,6))/out.T;
% out.IPF = sin(atan(out.iq./out.id)-atan(out.fq./out.fd));
out.SOL = SOL;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uncomment next line to visualize motor geometries during
% optimization (slower)
% plotfemmproblem(FemmProblem),drawnow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cost functions
cost1 = -out.T;
cost2 = out.dTpu*100;
cost  = [cost1 cost2];

% delta = atan2(out.fq,out.fd) * 180/pi;
% geo.power_factor = mean(cosd(90+delta-gamma));

% penalize the solutions which are out of the expected range of interest
% max_exp_ripple, min_exp_torque
if strcmp(eval_type,'MO_OA')||strcmp(eval_type,'MO_GA')
    if (cost(2)>per.max_exp_ripple || cost(1)>-per.min_exp_torque)
        cost(1)=cost(1)*0.1;
        cost(2)=cost(2)*10;
    end
end

if strcmp(geo.RemoveTMPfile,'ON')
    cd ..
    rmdir(dirName,'s')
end
cd(currentDir)

