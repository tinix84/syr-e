function DrawMachineScript(dataSet,pathname,filename)

if nargin()<3
    pathname = [cd '\'];
    filename = 'newmachine.fem';
    disp('You do not insert pathname and filename.')
    disp('Default parameters:')
    disp(['pathname : ' pathname]);
    disp(['filename : ' filename]);
    disp('dataSet inserted from manual_dataSet.mat');
    
    manual_dataSet;
end


if ~strcmp(filename(end-3),'.')
    filename = [filename '.mat'];
end
dataSet.currentfilename = [filename(1:end-4) '.mat'];
if ~isequal(pathname(end),'\')
    pathname=[pathname '\'];
end
dataSet.currentpathname = pathname;



%% from DrawPushMachine.m

[bounds, ~, geo, per, mat] = data0(dataSet);
RQ = buildDefaultRQ(bounds);

% RQ defines the candidate machine 
% geo.pathname = cd;

[geo,~,mat] = interpretRQ(RQ,geo,mat);

% FemmProblem = loadfemmfile([currentDir filesep 'empty_case.fem']);
FemmProblem.ProbInfo.Frequency = 0;
FemmProblem.ProbInfo.Precision = 1e-8;
FemmProblem.ProbInfo.MinAngle = 15;
FemmProblem.ProbInfo.LengthUnits = 'millimeters';
FemmProblem.ProbInfo.Depth = geo.l;
FemmProblem.Segments = [];
FemmProblem.ArcSegments = [];
FemmProblem.Nodes = [];
FemmProblem.BoundaryProps = [];
FemmProblem.Circuits = [];
FemmProblem.BlockLabels = [];
FemmProblem.PointProps = [];
eval_type = 'singt';
openfemm
[geo,mat] = draw_motor_in_FEMM(geo,eval_type,mat);
delete('mot0.fem');

%% SAVE THE FILE ==========================================================
% writefemmfile(filename, FemmProblem);
mi_saveas([pathname filename(1:end-4) '.fem']);
mi_close, closefemm
geo.RQ = RQ;
dataSet.RQ = RQ;
% dataSet.RQnames = geo.RQnames;
% dataSet.Dalpha1BouCheck = 0;
% dataSet.DalphaBouCheck = 0;
% dataSet.hcBouCheck = 0;
% dataSet.DxBouCheck = 0;
% dataSet.GammaBouCheck = 0;
% dataSet.GapBouCheck  = 0;
% dataSet.BrBouCheck  = 0;
% dataSet.AirgapRadiusBouCheck  = 0;
% dataSet.ToothWidthBouCheck  = 0;
% dataSet.ToothLengthBouCheck  = 0;
% dataSet.StatorSlotOpenBouCheck  = 0;
% dataSet.ToothTangDepthBouCheck  = 0;

filename = strrep(filename,'fem','mat');
save([pathname filename],'FemmProblem','geo','per','dataSet','mat');

disp(' ')
disp([filename ' saved in ' pathname])



