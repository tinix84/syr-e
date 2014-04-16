
%% Overall Description
% This code implements a basic multi-objective optimization algorithm based
% on Diferential Evolution (DE) algorithm:
%
% Storn, R., Price, K., 1997. Differential evolution: A simple and 
% efficient heuristic for global optimization over continuous spaces. 
% Journal of Global Optimization 11, 341 – 359.
%
% When one objective is optimized, the standard DE runs; if two or more
% objectives are optimized, the greedy selection step in DE algorithm is 
% performed using a dominance relation.
%%

close all;
clc;
data0
%% Variables regarding the optimization problem

dat.FieldD =bounds; % Initialization bounds
dat.Initial=bounds; % Optimization bounds (see data0.m)
dat.NOBJ = 2;                          % Number of objectives
dat.NRES = 0;                          % Number of constraints
dat.NVAR   = size(bounds,1);                       % Numero of decision variables
dat.mop = str2func('evaluateF');        % Cost function
%dat.CostProblem=@ZDT1;                 % Cost function instance

%%%%%%%%%% FEMM fitness handle %%%%%%%%%%%%%%%%%%%%%%%%%%
eval_type='MO_OA';
FitnessFunction = @(x)FEMMfitness(x,eval_type);
dat.CostProblem = FitnessFunction;                 % Cost function instance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Variables regarding the optimization algorithm
% For guidelines for the parameter tuning see:
%
% Storn, R., Price, K., 1997. Differential evolution: A simple and 
% efficient heuristic for global optimization over continuous spaces. 
% Journal of Global Optimization 11, 341 – 359.
%
% Das, S., Suganthan, P. N., 2010. Differential evolution: A survey of the 
% state-of-the-art. IEEE Transactions on Evolutionary Computation. Vol 15,
% 4 - 31.
%
dat.XPOP = 36;             % Population size
dat.Esc = 0.5;                         % Scaling factor
dat.Pm= 0.2;                           % Croosover Probability
%
%% Other variables
%
dat.InitialPop=[];                     % Initial population (if any)
dat.MAXGEN =300;                     % Generation bound
dat.MAXFUNEVALS = 150*dat.NVAR...  % Function evaluations bound
    *dat.NOBJ;                         
dat.SaveResults='yes';                 % Write 'yes' if you want to 
                                           % save your results after the
                                           % optimization process;
                                           % otherwise, write 'no';
%% Initialization (don't modify)
dat.CounterGEN=0;
dat.CounterFES=0;
%% Put here the variables required by your code (if any).
%
%
%
%% 
%
OUT=MODE(dat);                         % Run the algorithm.
%
%% Release and bug report:
%
% November 2012: Initial release
