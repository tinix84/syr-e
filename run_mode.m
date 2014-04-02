
%% Main Script %%
% runs the Multi-Objective optimization algorithm based on Differential
% Evolution

close all; clc;
format long
% FEMM must be installed in the defalut directory c:\
% FEMM can be downloaded at www.femm.info

warning('off','MATLAB:dispatcher:InexactMatch')

% load the file containing the input data 
data0;

eval_type = 'MO_OA';

FitnessFunction = @(x)FEMMfitness(x,eval_type);

% save geo_mot0 geo per mat
% load geo_mot_temp;
%%%%%%%%%%%%%%

if matlabpool('Size')>0
    matlabpool close force
end
matlabpool local
% matlabpool(5)

[x,fval,Pareto_front, Pareto_Fvals,exitFlag,output] = ...
    mode_optim(FitnessFunction,5,bounds(:,1),bounds(:,2),...
    'MaxFunEvals',100,'display','on','dispPLOT','on');

%% save the results
case_name = ['end_' datestr(now)];
case_name(case_name == ' ') = '_';
case_name(case_name == ':') = '-';
save(['results\' case_name]);

matlabpool close