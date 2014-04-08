%% MODE
% Multi-objective Evolutionary Algorithm (MOEA) based on Differential
% Evolution (DE).
% It implements a greedy selection based on pure dominance.
% DE algorithm has been introduced in:
%
% Storn, R., Price, K., 1997. Differential evolution: A simple and 
% efficient heuristic for global optimization over continuous spaces. 
% Journal of Global Optimization 11, 341 – 359.
%%


function OUT=MODE(options)

%% Reading parameters from options
generations     = options.MAXGEN;    % Maximum number of generations.
Xpop            = options.XPOP;      % Population size.
Nvar            = options.NVAR;      % Number of decision variables.
Nobj            = options.NOBJ;      % Number of objectives.
Bounds          = options.FieldD;    % Optimization bounds.
Initial         = options.Initial;   % Initialization bounds.
scaleFactor     = options.Esc;       % Scaling fator in DE algorithm.
crossOverRate   = options.Pm;        % Crossover probability in DE algorithm.
mop             = options.mop;       % Cost function.



%% Initial random population
Parent = zeros(Xpop,Nvar);  % Parent population.
Mutant = zeros(Xpop,Nvar);  % Mutant population.
Child  = zeros(Xpop,Nvar);  % Child population.
FES    = 0;                 % Function Evaluation.

for xpop=1:Xpop
    for nvar=1:Nvar
        Parent(xpop,nvar) = Initial(nvar,1)+(Initial(nvar,2) - Initial(nvar,1))*rand();
    end
end

if size(options.InitialPop,1)>=1
    Parent(1:size(options.InitialPop,1),:)=options.InitialPop;
end
disp('------------------------------------------------')
disp(['Initialization process']);
disp(['Evaluating initial population']);
disp('------------------------------------------------')

JxParent = mop(Parent,options);
FES = FES+Xpop;   

%% Evolution process

for n=1:generations 
    
    for xpop=1:Xpop
        rev=randperm(Xpop);
        
        %% Mutant vector calculation
        Mutant(xpop,:)= Parent(rev(1,1),:)+scaleFactor*(Parent(rev(1,2),:)-Parent(rev(1,3),:));
        
        for nvar=1:Nvar %Bounds are always verified
            if Mutant(xpop,nvar)<Bounds(nvar,1)
                Mutant(xpop,nvar) = Bounds(nvar,1);
            elseif Mutant(xpop,nvar)>Bounds(nvar,2)
                Mutant(xpop,nvar)=Bounds(nvar,1);
            end
        end
        
        %% Crossover operator
        for nvar=1:Nvar
            if rand() > crossOverRate
                Child(xpop,nvar) = Parent(xpop,nvar);
            else
                Child(xpop,nvar) = Mutant(xpop,nvar);
            end
        end

    end

    JxChild = mop(Child,options);
    FES=FES+Xpop;

    %% Selection
    for xpop=1:Xpop
        if JxChild(xpop,:) <= JxParent(xpop,:) 
            Parent(xpop,:) = Child(xpop,:);
            JxParent(xpop,:) = JxChild(xpop,:);
        end
    end
    
	PFront=JxParent;
	PSet=Parent;

    OUT.Xpop           = Parent;   % Population
    OUT.Jpop           = JxParent; % Poopulation's Objective Vector
    OUT.PSet           = PSet;     % Pareto Set
    OUT.PFront         = PFront;   % Pareto Front
    OUT.Param          = options;  % MODE Parameters
    options.CounterGEN = n;
    options.CounterFES = FES;
    
    [OUT, options]=PrinterDisplay(OUT,options); % To print results on screen
    
    if FES>options.MAXFUNEVALS || n>options.MAXGEN
        disp('Termination criteria reached.')
        break;
    end
end

OUT.Xpop=PSet;
OUT.Jpop=PFront;
[OUT.PFront, OUT.PSet]=DominanceFilter(PFront,PSet); %A Dominance Filter


if strcmp(options.SaveResults,'yes')
    save(['OUT_' datestr(now,30)],'OUT'); %Results are saved
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
for xpop=1:size(F,1)
    if Nobj==1
         hold on;
        plot(options.CounterGEN,log(min(F(:,1))),'dk', ...
            'MarkerFaceColor','k'); grid on; hold on;
    elseif Nobj==2
        hold on;
        plot(F(xpop,1),F(xpop,2),'dk','MarkerFaceColor','k');...
            grid on; hold on;
    elseif Nobj==3
         hold on;
        plot3(F(xpop,1),F(xpop,2),F(xpop,3),'dk','MarkerFaceColor','k');...
            grid on; hold on;
    end
end

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

%% Dominance Filter
%
% A filter based on dominance criteria
%
function [F, C]=DominanceFilter(F,C)

Xpop=size(F,1);
Nobj=size(F,2);
Nvar=size(C,2);
k=0;

for xpop=1:Xpop
    dom=0;
    
    for i=1:Xpop
        if F(xpop,:)==F(i,:)
            if xpop > i
                dom=1;
                break;
            end
        else
            if F(xpop,:)>=F(i,:)
                dom=1;
                break;
            end
        end
    end
    
    if dom==0
        k=k+1;
        F(k,:)=F(xpop,:);
        C(k,:)=C(xpop,:);
    end
end
F=F(1:k,:);
C=C(1:k,:);

%% Release and bug report:
%
% November 2012: Initial release