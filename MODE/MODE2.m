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

%% MODE
% Multi-objective Evolutionary Algorithm (MOEA) based on Differential
% Evolution (DE).
% It implements a greedy selection based on pure dominance.
% DE algorithm has been introduced in:
%
% Storn, R., Price, K., 1997. Differential evolution: A simple and
% efficient heuristic for global optimization over continuous spaces.
% Journal of Global Optimization 11, 341 � 359.
%%


function OUT=MODE2(options, dataSet)

%% Reading parameters from options
generations     = options.MAXGEN;    % Maximum number of generations.
populationSize  = options.XPOP;      % Population size.
numVariables    = options.NVAR;      % Number of decision variables.
numObjectives   = options.NOBJ;      % Number of objectives.
Bounds          = options.FieldD;    % Optimization bounds.
Initial         = options.Initial;   % Initialization bounds.
scaleFactor     = options.Esc;       % Scaling fator in DE algorithm.
crossOverRate   = options.Pm;        % Crossover probability in DE algorithm.
costFunction    = options.mop;       % Cost function.
OUT.MatrixPop   = [];
OUT.MatrixFitness = [];
OUT.MatrixPFront  = [];
OUT.MatrixPset    = [];
tau1            = options.tau1;

%% Initial random population
Parent = zeros(populationSize,numVariables);  % Parent population.
Mutant = zeros(populationSize,numVariables);  % Mutant population.
Child  = zeros(populationSize,numVariables);  % Child population.
functionEvaluations    = 0;                 % Function Evaluation.

for i=1:populationSize
    for nvar=1:numVariables
        Parent(i,nvar) = Initial(nvar,1)+(Initial(nvar,2) - Initial(nvar,1))*rand();
    end
end

if size(options.InitialPop,1)>=1
    Parent(1:size(options.InitialPop,1),:)=options.InitialPop;
end
disp('------------------------------------------------')
disp(['Initialization process']); %#ok<*NBRAK>
disp(['Evaluating initial population']);
disp('------------------------------------------------')

JxParent = costFunction(Parent,options);
functionEvaluations = functionEvaluations+populationSize;

%% Evolution process

for n=1:generations
    for i=1:populationSize
        rev=randperm(populationSize);
        
        %% Mutant vector calculation
        if rand<tau1
            sf=rand;
        else
            sf=scaleFactor;
        end
        Mutant(i,:)= Parent(rev(1),:)+sf*(Parent(rev(2),:)-Parent(rev(3),:));
        
        for nvar=1:numVariables %Bounds are always verified
            if Mutant(i,nvar)<Bounds(nvar,1)
                Mutant(i,nvar) = Bounds(nvar,1);
            elseif Mutant(i,nvar)>Bounds(nvar,2)
                Mutant(i,nvar)=Bounds(nvar,1);
            end
        end
        
        %% Crossover operator
        for nvar=1:numVariables
            if rand > crossOverRate
                Child(i,nvar) = Parent(i,nvar);
            else
                Child(i,nvar) = Mutant(i,nvar);
            end
        end
        
    end
    disp('------------------------------------------------')
    disp(['Evaluating ' num2str(n) '-th generation...']);
    disp('------------------------------------------------')
    
    JxChild = costFunction(Child,options);
    functionEvaluations=functionEvaluations+populationSize;
    
    %% Selection
    for k=1:populationSize
        if JxChild(k,:) <= JxParent(k,:)
            Parent(k,:) = Child(k,:);
            JxParent(k,:) = JxChild(k,:);
        end
    end
    tempPopulation = [Parent JxParent; Child JxChild];
    tempPopulation = nonDominationSort(tempPopulation, numObjectives, numVariables);
    tempPopulation = cutPopulation(tempPopulation, numObjectives, numVariables,populationSize);
    
    Parent     = tempPopulation(:,1:numVariables);
    Objectives = tempPopulation(:,numVariables+1:numVariables+numObjectives);
    Ranks      = tempPopulation(:,numVariables+1+numObjectives);
    %CrowdDist  = tempPopulation(:,end);
    
    PFront = Objectives(Ranks==1,:);
    PSet   = Parent(Ranks==1,:);
    OUT.MatrixPop      = cat(3, OUT.MatrixPop, Parent);
    OUT.MatrixFitness  = cat(3, OUT.MatrixFitness, Objectives);
    OUT.Xpop           = Parent;   % Population
    OUT.Jpop           = Objectives; % Poopulation's Objective Vector
    OUT.PSet           = PSet;     % Pareto Set
    OUT.PFront         = PFront;   % Pareto Front
    OUT.MatrixPset{n}  = PSet;
    OUT.MatrixPFront{n}   = PFront;
    OUT.Param          = options;  % MODE Parameters
    OUT.eval_type      = options.eval_type;
    options.CounterGEN = n;
    options.CounterFES = functionEvaluations;
    
    [OUT, options]=PrinterDisplay(OUT,options); % To print results on screen
    
    if functionEvaluations>options.MAXFUNEVALS || n>options.MAXGEN
        disp('Termination criteria reached.')
        break;
    end
end

OUT.Xpop=Parent;
OUT.Jpop=Objectives;
%[OUT.PFront, OUT.PSet]=DominanceFilter(PFront,PSet); %A Dominance Filter
s=size(OUT.MatrixPop);

M1=[];
M2=[];
if numel(s)<3
    M1=[M1;OUT.MatrixPset{1}]; % #ok<AGROW>
    M2=[M2;OUT.MatrixPFront{1}]; % #ok<AGROW>
else
    for i=1:s(3)
        M1=[M1;OUT.MatrixPset{i}]; %#ok<AGROW>
        M2=[M2;OUT.MatrixPFront{i}]; %#ok<AGROW>
    end
end
I=isparetosetMember(M2);
OUT.PFront=(M2(I==1,:));
OUT.PSet=(M1(I==1,:));

if strcmp(options.SaveResults,'yes')
    geo0=options.geo0;
    per=options.per;
    mat=options.mat;
    thisfilepath = fileparts(which('data0.m'));
    filename=fullfile(thisfilepath,'results',['OUT_' datestr(now,30)]);
    save(filename,'OUT','per','geo0','dataSet','mat'); %Results are saved
    clear geo0 per
end

disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp('Red  asterisks : Set Calculated.')
disp('Black diamonds : Filtered Set.')
if strcmp(options.SaveResults,'yes')
    disp(['Check out OUT_' datestr(now,30) ...
        ' variable on folder for results.'])
end
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')

F=OUT.PFront;

plot(F(:,1),F(:,2),'dk','MarkerFaceColor','k');...
    grid on; hold on;
% PostProcessing of current optimization result
% postProcExecution = questdlg('Do you want to run post processing??','Postprocessing','Yes','No','Yes');
% if strcmp(postProcExecution,'Yes')
disp('PostProcessing of current optimization result...');
[~,ff]=fileparts(filename);
if strcmp(func2str(options.CostProblem), '@(x)FEMMfitness(x,geo,per,mat,eval_type)')
    evalParetoFront(ff,dataSet);
else
    evalParetoFrontX(ff,dataSet);
end

% end

%% Print and Display information
% Modify at your convenience
%
function [OUT, Dat]=PrinterDisplay(OUT,Dat)

disp('------------------------------------------------')
disp(['Generation: ' num2str(Dat.CounterGEN)]);
disp(['FEs: ' num2str(Dat.CounterFES)]);
disp(['Pareto Front Size: ' mat2str(size(OUT.PFront,1))]);
disp('------------------------------------------------')

if mod(Dat.CounterGEN,1)==0
    if Dat.NOBJ==3
        
        plot3(OUT.PFront(:,1),OUT.PFront(:,2),OUT.PFront(:,3),'*r');
        grid on;hold on,drawnow
    elseif Dat.NOBJ==2
        
        plot(OUT.PFront(:,1),OUT.PFront(:,2),'*r'); grid on;hold on,drawnow
    elseif Dat.NOBJ==1
        
        plot(Dat.CounterGEN,log(min(OUT.PFront(:,1))),'*r'); ...
            grid on; hold on,drawnow
    end
end


%% Release and bug report:
%
% November 2012: Initial release