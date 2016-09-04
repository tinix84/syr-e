function [dataSet,geo,per,mat] = back_compatibility(dataSet,geo,per)

flag = 0;

%% from version 1.1 to version 1.4 (november 2015)
% disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')
% disp('update to version 1.4')
if ~isfield(dataSet,'TargetCopperTemp')
    dataSet.TargetCopperTemp = dataSet.CopperTemp;
    dataSet=rmfield(dataSet,'CopperTemp');
    disp('v1.4 - renamed CopperTemp')
    flag = 1;
end

if ~isfield(dataSet,'HousingTemp')
    dataSet.HousingTemp = 50;
    disp('v1.4 - added housing temperature')
    flag = 1;
end

if ~isfield(dataSet,'EstimatedCopperTemp')
    dataSet.EstimatedCopperTemp = per.tempcu;
    disp('v1.4 - added estimated copper temperature')
    flag = 1;
end

if ~isfield(dataSet,'randFactor')
    dataSet.randFactor = 0;
    flag = 1;
end

if ~isfield(dataSet,'MagLoadingYoke')
    dataSet.MagLoadingYoke = 0.5;
    disp('v1.4 - added magnetic loading of the yoke')
    flag = 1;
end

if ~isfield(dataSet,'MagLoadingTooth')
    dataSet.MagLoadingTooth = 1;
    disp('v1.4 - added magnetic loading of the tooth')
    flag = 1;
end

if ~isfield(dataSet,'ThicknessOfPM')
    dataSet.ThicknessOfPM = geo.g*6;
    disp('v1.4 - added thickness of PM (SPM)')
    flag = 1;
end

if ~isfield(dataSet,'AngleSpanOfPM')
    dataSet.AngleSpanOfPM = 150;
    disp('v1.4 - added angle span of PM (SPM)')
    flag = 1;
end

if ~isfield(dataSet,'DepthOfBarrier')
    dataSet.DepthOfBarrier = ones(1,geo.nlay);
    disp('v1.4 - added depth of barrier')
    flag = 1;
end

if ~isfield(dataSet,'StatorSlotOpenBou')
    dataSet.StatorSlotOpenBou = [0.2 0.3];
    dataSet.StatorSlotOpenBouCheck = 0;
    disp('v1.4 - added boundaries for acs optimization')
    flag = 1;
    
end

if ~isfield(dataSet,'ToothTangDepthBou')
    dataSet.ToothTangDepthBou = geo.g*[1 3];
    dataSet.ToothTangDepthBouCheck = 0;
    disp('v1.4 - added boundaries for ttd optimization')
    flag = 1;
end

%% check the dimension of rotor parameters

if length(dataSet.HCpu)~=dataSet.NumOfLayers
    dataSet.HCpu = dataSet.HCpu(1:dataSet.NumOfLayers);
    disp('v1.4 - correct dataSet.HCpu')
    flag = 1;
end

if length(dataSet.ALPHApu)~=dataSet.NumOfLayers
    dataSet.ALPHApu = dataSet.ALPHApu(1:dataSet.NumOfLayers);
    disp('v1.4 - correct dataSet.ALPHApu')
    flag = 1;
end

% if length(dataSet.DepthOfBarrier)~=dataSet.NumOfLayers
%     dataSet.DepthOfBarrier = geo.dx;
%     dataSet.DepthOfBarrier = dataSet.DepthOfBarrier(1:dataSet.NumOfLayers);
%     disp('v1.4 - correct dataSet.DepthOfBarrier')
%     flag = 1;
% end

%% check dimension of RQ

if length(dataSet.RQnames)~=length(dataSet.RQ)
    dataSet.RQ = geo.RQ;
    disp('v1.4 - correct dataSet.RQ')
    flag = 1;
end


%% from version 1.4 to rev 260
% disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')
% disp('update to rev 260')

if ~isfield(dataSet,'LossEvaluationCheck')
    dataSet.LossEvaluationCheck = 0;
    dataSet.HysteresisLossFactor = 0;
    dataSet.HysteresisFrequencyFactor = 0;
    dataSet.HysteresisFluxDenFactor = 0;
    dataSet.EddyCurLossFactor = 0;
    dataSet.EddyCurLossFactorEdit = 0;
    dataSet.EvalSpeed = 0;
    dataSet.IronMassDen = geo.rhoFE;
    disp('rev260 - added parameters for iron losses evaluation')
    flag = 1;
end

if ~isfield(dataSet,'TorqueOptCheck')
    dataSet.TorqueOptCheck = 1;
    dataSet.TorRipOptCheck = 1;
    disp('rev260 - added check for optimization objective')
    flag = 1;
end

if ~isfield(dataSet,'Qs')
    Q = 6*geo.p*geo.q;
    t2 = gcd(round(dataSet.NumOfSlots*6*dataSet.NumOfPolePairs),2*dataSet.NumOfPolePairs);
    dataSet.Qs = Q/t2;
    clear t2;
    disp('rev260 - added number of simulated stator slot')
    flag = 1;
end

if ~isfield(dataSet,'xRange')
    dataSet.xRange = [0.5 0.8];
    dataSet.bRange = [0.3 0.7];
    disp('rev260 - added (x,b) range for syrmDesign')
    flag = 1;
end

if ~isfield(dataSet,'BarFillFac')
    dataSet.BarFillFac = 0;
    disp('rev260 - added barrier filling factor')
    flag = 1;
end



%% from rev 260 to 261
% disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')
% disp('update to 261')

if ~isfield(dataSet,'MassCuOptCheck')
    dataSet.MassCuOptCheck = 0;
    dataSet.MaxCuMass = 0;
    disp('rev261 - added copper mass in optimization objective')
    flag = 1;
end

if ~isfield(dataSet,'SlotLayerPosCheck')
    dataSet.SlotLayerPosCheck = 0;
    disp('rev261 - added slot layer position in GUI')
    flag = 1;
end

if ~isfield(dataSet,'RadRibCheck')
    dataSet.RadRibCheck = 0;
    dataSet.RadRibEdit = zeros(1,geo.nlay);
    disp('rev261 - added radial ribs in GUI')
    flag = 1;
end

%% rewriting geo and per
[~, ~, geo, per, mat] = data0(dataSet);

%% message in command window if some data are added
if flag
%     disp('The selected project is in an older version of SyR-e: save machine to update to the latest version')
    disp('...');
end