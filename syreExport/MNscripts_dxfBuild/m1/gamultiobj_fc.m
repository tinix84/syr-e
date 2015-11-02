function [x,fval,exitFlag,output,population,scores] = gamultiobj(fun,nvars,Aineq,bineq,Aeq,beq,lb,ub,options)
%GAMULTIOBJ Multiobjective optimization using genetic algorithm.
%   GAMULTIOBJ attempts to solve multiobjective problems of the form:
%       min F(X)  subject to:  A*X <= b, Aeq*X = beq (linear constraints)
%        X                     lb <= X <= ub (bound constraints)
%
%   X = GAMULTIOBJ(FITNESSFCN,NVARS) finds a local Pareto set X of the
%   objective functions defined in FITNESSFCN. NVARS is the dimension of
%   the optimization problem (number of decision variables). FITNESSFCN
%   accepts a vector X of size 1-by-NVARS and returns a vector of
%   size 1-by-numberOfObjectives evaluated at a decision variable. X is
%   a matrix with NVARS columns. The number of rows in X is the same as the
%   number of Pareto solutions. All solutions in a Pareto set are equally
%   optimal, and it is up to the designer to select a solution in the Pareto
%   set depending on the application.
%
%   X = GAMULTIOBJ(FITNESSFCN,NVARS,A,b) finds a local Pareto set X of the
%   objective functions defined in FITNESSFCN, subject to the linear
%   inequalities A*X <= B. Linear constraints are supported only for
%   default PopulationType option ('doubleVector'); other population types
%   e.g., 'bitString' and 'custom' are not supported.
%
%   X = GAMULTIOBJ(FITNESSFCN,NVARS,A,b,Aeq,beq) finds a local Pareto set X
%   of the objective functions defined in FITNESSFCN, subject to the linear
%   equalities Aeq*X = beq as well as the linear inequalities A*X <= b. (Set
%   A=[] and b=[] if no inequalities exist.) Linear constraints are supported
%   only for default PopulationType option ('doubleVector'); other population
%   types e.g., 'bitString' and 'custom' are not supported.
%
%   X = GAMULTIOBJ(FITNESSFCN,NVARS,A,b,Aeq,beq,lb,ub) defines a set of
%   lower and upper bounds on the design variables, X, so that a local Pareto
%   set is found in the range lb <= X <= ub. Use empty matrices for lb and ub
%   if no bounds exist. Set lb(i) = -Inf if X(i) is unbounded below;  set
%   ub(i) = Inf if X(i) is unbounded above. Bound constraints are
%   supported only for default PopulationType option ('doubleVector');
%   other population types e.g., 'bitString' and 'custom' are not
%   supported.
%
%   X = GAMULTIOBJ(FITNESSFCN,NVARS,A,b,Aeq,beq,lb,ub,options) finds a
%   Pareto set X with the default optimization parameters replaced by
%   values in the structure OPTIONS. OPTIONS can be created with the
%   GAOPTIMSET function. See GAOPTIMSET for details.
%
%   X = GAMULTIOBJ(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure that has the following fields:
%       fitnessfcn: <Fitness function>
%            nvars: <Number of design variables>
%            Aineq: <A matrix for inequality constraints>
%            bineq: <b vector for inequality constraints>
%              Aeq: <Aeq matrix for equality constraints>
%              beq: <beq vector for equality constraints>
%               lb: <Lower bound on X>
%               ub: <Upper bound on X>
%          options: <Options structure created with GAOPTIMSET>
%           solver: <solver name 'gamultiobj'>
%         rngstate: <State of the random number generator>
%
%   [X,FVAL] = GAMULTIOBJ(FITNESSFCN,NVARS, ...) in addition returns a
%   matrix
%   FVAL, the value of all the objective functions defined in FITNESSFCN at
%   all the solutions in X. FVAL has numberOfObjectives columns and same number
%   of rows as does X.
%
%   [X,FVAL,EXITFLAG] = GAMULTIOBJ(FITNESSFCN,NVARS, ...) in addition returns
%   EXITFLAG which describes the exit condition of GAMULTIOBJ. Possible values
%   of EXITFLAG and the corresponding exit conditions are
%
%     1 Average change in value of the spread of Pareto set over
%        options.StallGenLimit generations less than options.TolFun.
%     0 Maximum number of generations exceeded.
%    -1 Optimization terminated by the output or plot function.
%    -2 No feasible point found.
%    -5 Time limit exceeded.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = GAMULTIOBJ(FITNESSFCN,NVARS, ...) in addition
%   returns a structure OUTPUT with the following fields:
%            rngstate: <State of the random number generator before GA started>
%         generations: <Total number of generations, excluding HybridFcn iterations>
%           funccount: <Total number of function evaluations>
%       maxconstraint: <Maximum constraint violation>, if any
%             message: <GAMULTIOBJ termination message>
%
%   [X,FVAL,EXITFLAG,OUTPUT,POPULATION] = GAMULTIOBJ(FITNESSFCN, ...) in
%   addition returns the final POPULATION at termination. 
%
%   [X,FVAL,EXITFLAG,OUTPUT,POPULATION,SCORE] = GAMULTIOBJ(FITNESSFCN, ...)
%   in addition returns the SCORE of the final POPULATION.
%
%
%   Example:
%    Multiobjective minimization of two functions 'ackleyfcn' and 'shufcn'
%
%      fun1 = @(x) ackleyfcn(x);
%      fun2 = @(x) shufcn(x);
%      % Combine two objectives 'fun1' and 'fun2' 
%      fun1and2 = @(x) [fun1(x) fun2(x)];
%      % Bound constraints on X
%      lb = [-10 -10]; ub = [10 10];
%      % Specify the initial range for population
%      options = gaoptimset('PopInitRange',[lb;ub]);
%      % Minimize using GAMULTIOBJ
%      [x,fval] = gamultiobj(fun1and2,2,[],[],[],[],lb,ub,options)
%
%    Display Pareto fronts and rank of individuals while GAMULTIOBJ
%    minimizes 
%
%      options = gaoptimset('PopInitRange',[lb;ub]);
%      options = gaoptimset(options,'PlotFcns',{@gaplotpareto,@gaplotrankhist});
%      [x,fval,exitflag,output] = gamultiobj(fun1and2,2,[],[],[],[],lb,ub,options)
%
%
%   See also GAOPTIMSET, PATTERNSEARCH, GA, @.

%   Reference: Kalyanmoy Deb, "Multi-Objective Optimization using
%   Evolutionary Algorithms", John Wiley & Sons ISBN 047187339


%   Copyright 2007-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $  $Date: 2008/05/23 15:33:39 $

% If the first arg is not a gaoptimset, then it's a fitness function followed by a genome
% length. Here we make a gaoptimset from the args.
defaultopt = struct('PopulationType', 'doubleVector', ...
    'PopInitRange', [0;1], ...
    'PopulationSize', '15*numberOfVariables', ...
    'CrossoverFraction', 0.8, ...
    'ParetoFraction', 0.35, ...
    'MigrationDirection','forward', ...
    'MigrationInterval',20, ...
    'MigrationFraction',0.2, ...
    'Generations', '200*numberOfVariables', ...
    'TimeLimit', inf, ...
    'StallGenLimit', 100, ...
    'TolFun', 1e-4, ...
    'TolCon', 1e-6, ...
    'InitialPopulation',[], ...
    'InitialScores', [], ...
    'PlotInterval',1, ...
    'CreationFcn',@gacreationuniform, ...
    'SelectionFcn', {{@selectiontournament,2}}, ...
    'CrossoverFcn',@crossoverintermediate, ...
    'MutationFcn',@mutationadaptfeasible, ...
    'DistanceMeasureFcn',{{@distancecrowding, 'phenotype'}}, ...
    'HybridFcn',[], ...
    'Display', 'final', ...
    'PlotFcns', [], ...
    'OutputFcns', [], ...
    'Vectorized', 'off', ...
    'UseParallel', 'never');

% If just 'defaults' passed in, return the default options in X
if nargin == 1 && nargout <= 1 && isequal(fun,'defaults')
    x = defaultopt;
    return
end

% Check number of input arguments
errmsg = nargchk(1,9,nargin);
if ~isempty(errmsg)
    error('gads:gamultiobj:numberOfInputs',[errmsg,' GAMULTIOBJ requires at least 1 input argument.']);
end

if nargin < 9,  options = [];
    if nargin < 8, ub = [];
        if nargin < 7, lb = [];
            if nargin <6, beq = [];
                if nargin <5, Aeq = [];
                    if nargin < 4, bineq = [];
                        if nargin < 3, Aineq = [];
                        end
                    end
                end
            end
        end
    end
end

% One input argument is for problem structure
if nargin == 1
    if isa(fun,'struct')
        [fun,nvars,Aineq,bineq,Aeq,beq,lb,ub,rngstate,options] = separateOptimStruct(fun);
        % Reset the random number generators
        resetDfltRng(rngstate);
    else % Single input and non-structure.
        error('gads:gamultiobj:invalidStructInput','The input should be a structure with valid fields or provide at least two arguments to GAMULTIOBJ.' );
    end
end

% If fun is a cell array with additional arguments get the function handle
if iscell(fun)
    FitnessFcn = fun{1};
else
    FitnessFcn = fun;
end
% Only function handles or inlines are allowed for FitnessFcn
if isempty(FitnessFcn) ||  ~(isa(FitnessFcn,'inline') || isa(FitnessFcn,'function_handle'))
    error('gads:gamultiobj:needFunctionHandle','Fitness function must be a function handle.');
end

% We need to check the nvars here before we call any solver
valid =  isnumeric(nvars) && isscalar(nvars)&& (nvars > 0) ...
    && (nvars == floor(nvars));
if ~valid
    error('gads:gamultiobj:notValidNvars','Number of variables (NVARS) must be a positive integer.');
end
user_options = options;
% Use default options if empty
if ~isempty(options) && ~isa(options,'struct')
        error('gads:gamultiobj:optionsNotAStruct','Ninth input argument must be a valid structure created with GAOPTIMSET.');
elseif isempty(options)
    options = defaultopt;
end
% Take defaults for parameters that are not in options structure
options = gaoptimset(defaultopt,options);

% All inputs should be double
try
    dataType = superiorfloat(nvars,Aineq,bineq,Aeq,beq,lb,ub);
    if ~isequal('double', dataType)
        error('gads:gamultiobj:dataType', ...
            'GAMULTIOBJ only accepts inputs of data type double.')
    end
catch
    error('gads:gamultiobj:dataType', ...
        'GAMULTIOBJ only accepts inputs of data type double.')
end

% Initialize unused variable for nonlcon
nonlcon = [];

% Validate constraints and options
[x,fval,exitFlag,output,population,scores,FitnessFcn,nvars,Aineq,bineq,Aeq,beq,lb,ub, ...
    NonconFcn,options] = gacommon(nvars,fun,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options,user_options);

if exitFlag < 0 % Infeasibility
    return;
end

% Call appropriate single objective optimization solver
[x,fval,exitFlag,output,population,scores] = gamultiobjsolve_fc(FitnessFcn,nvars, ...
     Aineq,bineq,Aeq,beq,lb,ub,options,output);

 
