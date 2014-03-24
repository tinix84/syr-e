
%% Main Script %%
% runs the Multi-Objective optimization algorithm based on Differential
% Evolution

clear all; close all; clc;

% FEMM must be installed in the defalut directory c:\
% FEMM can be downloaded at www.femm.info
addpath('d:\femm42_beta\mfiles\');
addpath('MODE');
% addpath('GODLIKE');
addpath mfiles;

warning('off','MATLAB:dispatcher:InexactMatch')

% load the file containing the input data 
data0;

global eval_type options_global

eval_type = 'MO_OA';

FitnessFunction = @FEMMfitness;

% save geo_mot0 geo per mat
% load geo_mot_temp;
%%%%%%%%%%%%%%

if matlabpool('Size')>0
    matlabpool close force
end
% matlabpool
% matlabpool(5)

[x,fval,Pareto_front, Pareto_Fvals,exitFlag,output] = ...
    mode_optim(FitnessFunction,5,bounds(:,1),bounds(:,2),...
    'MaxFunEvals',100,'display','on','dispPLOT','on');

%% save the results
case_name = ['end_' datestr(now)];
case_name(case_name == ' ') = '_';
case_name(case_name == ':') = '-';
save(['results\' case_name]);

