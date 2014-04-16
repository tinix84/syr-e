function [cost] = FEMMfitnessX(RQ,eval_type)

RQ=RQ';

geo.pathname=cd;
data0;                            

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

p = geo.p;                      % paia poli
nlay = geo.nlay;                % numero delle barriere

% transform the 
% note: the vector of the inputs RQ contains
% RQ(1) = dalpha(1) [mec deg]
% RQ(2) = dalpha(2) [p.u.]
% ...
% RQ(nlay) = dalpha(nlay); [p.u.]


% RQ(end) = gamma

first_index = 2;
last_index = first_index + nlay - 1;

dalpha_pu = RQ(first_index:last_index);

% if the sum of the pu angles is too large, it is scaled down     
if sum(dalpha_pu) > 1
    dalpha_pu = dalpha_pu/sum(dalpha_pu);
end

% dalpha(2) to dalpha(nlay) in degrees 
dalpha_temp = dalpha_pu * (90/p - RQ(1));
% all dalpha in mec degrees
geo.dalpha = [RQ(1) dalpha_temp(1:end-1)];

% SPESSORE DELLE BARRIERE: 'hc_pu'
first_index = last_index + 1;
last_index = first_index + nlay - 1;
size(RQ);
getComputerName();
geo.hc_pu = RQ(first_index:last_index);
% convert dalpha and hc_pu to alpha and hc
geo = calc_alpha_hc_delta_x0_2(geo);

if (strcmp(geo.RotType,'3U'))
    first_index = last_index + 1;
    last_index=first_index;
    geo.Delta_X=RQ(first_index:last_index);
elseif (strcmp(geo.RotType,'Fluid'))
    first_index = last_index + 1;
    last_index = first_index + nlay - 1;
    geo.Dfe=RQ(first_index:last_index);
end
% current phase angle
gamma = RQ(end);
mang=num2str(rand);
dirName=mang(3:end);
%while(~exist('tmp','dir'))
warning off MATLAB:MKDIR:DirectoryExists
    mkdir('tmp');
%end
cd('tmp')
while(exist(dirName,'dir'))
    mang=num2str(rand);
    dirName=mang(3:end);
end
mkdir(dirName);
cd(dirName);
copyfile('c:\empty_case.fem','.');
%openfemm
FemmProblem = newproblem_mfemm('planar', ...
                               'Frequency', 0, ...
                               'LengthUnits', 'millimeters',...
                               'Depth',geo.l,'MinAngle',25);
draw_motor_in_XFEMM
save geo_mot_temp      
% current amplitude used for the simulations
io = per.overload*calc_io(geo,per);
geo.io=io;
% current value use for FEMM simulation, only 1 turns is set in femm
io_femm=io*geo.Nbob;
% evaluates the candidate machine (T, DT, fd, fq)
[out,FemmProblem] = eval_motor_in_FEMMX(FemmProblem,geo,io_femm,gamma,eval_type,boundnameAPmove);

% Tn = out.Tn;
numsim = size(out.SOL,1);
% ripple_pu = out.ripple_pu;

% Cost functions to be minimized
cost1 = -out.Tn;
cost2 = out.ripple_pu*100;
cost  = [cost1 cost2];

delta = atan2(out.fq,out.fd) * 180/pi;
geo.power_factor = mean(cosd(delta-gamma));

% penalize the solutions which are out of the expected range of interest
% max_exp_ripple, min_exp_torque
if strcmp(eval_type,'MO_OA')
    if (cost(2)>per.max_exp_ripple || cost(1)>-per.min_exp_torque)
        cost(1)=cost(1)*0.1;
        cost(2)=cost(2)*10;
    end
end

save geo_mot_temp geo fem cost out RQ BLKLABELS per rotore2 statore

%closefemm
cd('..\..')

