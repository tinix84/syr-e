function [dataSet,geo,per,mat] = back_compatibility(dataSet,geo,per,Dflag)
% 
% [dataSet,geo,per,mat] = back_compatibility(dataSet,geo,per,Dflag)
% 
% Update the project to newer version of SyR-e
% INPUT : dataSet
%         geo
%         per
%         Dflag =1-->disp the modification / =0-->don't plot anything
% OUTPUT: dataSet
%         geo (geometry)
%         per (performance)
%         mat (material)

flag = 0;

if nargin()<4
    Dflag=1;
end

%% from version 1.1 to version 1.4 (november 2015)
% disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')
% disp('update to version 1.4')
if ~isfield(dataSet,'TargetCopperTemp')
    dataSet.TargetCopperTemp = dataSet.CopperTemp;
    dataSet=rmfield(dataSet,'CopperTemp');
    if Dflag
        disp('v1.4 - renamed CopperTemp')
    end
    flag = 1;
end

if ~isfield(dataSet,'HousingTemp')
    dataSet.HousingTemp = 50;
    if Dflag
        disp('v1.4 - added housing temperature')
    end
    flag = 1;
end

if ~isfield(dataSet,'EstimatedCopperTemp')
    dataSet.EstimatedCopperTemp = per.tempcu;
    if Dflag
        disp('v1.4 - added estimated copper temperature')
    end
    flag = 1;
end

if ~isfield(dataSet,'randFactor')
    dataSet.randFactor = 0;
    flag = 1;
end

if ~isfield(dataSet,'MagLoadingYoke')
    dataSet.MagLoadingYoke = 0.5;
    if Dflag
        disp('v1.4 - added magnetic loading of the yoke')
    end
    flag = 1;
end

if ~isfield(dataSet,'MagLoadingTooth')
    dataSet.MagLoadingTooth = 1;
    if Dflag
        disp('v1.4 - added magnetic loading of the tooth')
    end
    flag = 1;
end

if ~isfield(dataSet,'ThicknessOfPM')
    dataSet.ThicknessOfPM = geo.g*6;
    if Dflag
        disp('v1.4 - added thickness of PM (SPM)')
    end
    flag = 1;
end

if ~isfield(dataSet,'AngleSpanOfPM')
    dataSet.AngleSpanOfPM = 150;
    if Dflag
        disp('v1.4 - added angle span of PM (SPM)')
    end
    flag = 1;
end

if ~isfield(dataSet,'DepthOfBarrier')
    dataSet.DepthOfBarrier = ones(1,geo.nlay);
    if Dflag
        disp('v1.4 - added depth of barrier')
    end
    flag = 1;
end

if ~isfield(dataSet,'StatorSlotOpenBou')
    dataSet.StatorSlotOpenBou = [0.2 0.3];
    dataSet.StatorSlotOpenBouCheck = 0;
    if Dflag
        disp('v1.4 - added boundaries for acs optimization')
    end
    flag = 1;
    
end

if ~isfield(dataSet,'ToothTangDepthBou')
    dataSet.ToothTangDepthBou = geo.g*[1 3];
    dataSet.ToothTangDepthBouCheck = 0;
    if Dflag
        disp('v1.4 - added boundaries for ttd optimization')
    end
    flag = 1;
end

%% check the dimension of rotor parameters

if length(dataSet.HCpu)~=dataSet.NumOfLayers && ~strcmp(dataSet.TypeOfRotor,'SPM')
    dataSet.HCpu = dataSet.HCpu(1:dataSet.NumOfLayers);
    if Dflag
        disp('v1.4 - correct dataSet.HCpu')
    end
    flag = 1;
end

if length(dataSet.ALPHApu)~=dataSet.NumOfLayers && ~strcmp(dataSet.TypeOfRotor,'SPM')
    dataSet.ALPHApu = dataSet.ALPHApu(1:dataSet.NumOfLayers);
    if Dflag
        disp('v1.4 - correct dataSet.ALPHApu')
    end
    flag = 1;
end

if length(dataSet.DepthOfBarrier)~=dataSet.NumOfLayers && ~strcmp(dataSet.TypeOfRotor,'SPM')
    dataSet.DepthOfBarrier = geo.dx;
    dataSet.DepthOfBarrier = dataSet.DepthOfBarrier(1:dataSet.NumOfLayers);
    if Dflag
        disp('v1.4 - correct dataSet.DepthOfBarrier')
    end
    flag = 1;
end

%% check dimension of RQ

% [~,c]=size(dataSet.RQnames);
% if ~c
%     dataSet.RQ = [];
% end

% if length(dataSet.RQnames)~=length(dataSet.RQ)
%     if isfield(geo,'RQ') && isfield(geo,'RQnames')
%         dataSet.RQ=geo.RQ;
%         dataSet.RQnames=geo.RQnames;
%         if Dflag
%             disp('correct RQ and RQnames')
%         end
%     else
%         dataSet.RQ=[];
%         dataSet.RQnames=[];
%         dataSet.Dalpha1BouCheck=0;
%         dataSet.DalphaBouCheck=0;
%         dataSet.hcBouCheck=0;
%         dataSet.DxBouCheck=0;
%         dataSet.GammaBouCheck=0;
%         dataSet.GapBouCheck=0;
%         dataSet.BrBouCheck=0;
%         dataSet.AirgapRadiusBouCheck=0;
%         dataSet.ToothWidthBouCHeck=0;
%         dataSet.ToothLengthBouCheck=0;
%         dataSet.StatorSlotOpenBouCheck=0;
%         dataSet.ToothTangDepthBouCheck=0;
%         if Dflag
%             disp('error in RQ and RQnames: reset optimization variables')
%         end
%     end
             
%     dataSet.RQ = geo.RQ;
%     if Dflag
%         disp('v1.4 - correct dataSet.RQ')
%     end
%     flag = 1;
% end


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
    if Dflag
        disp('rev260 - added parameters for iron losses evaluation')
    end
    flag = 1;
end

if ~isfield(dataSet,'TorqueOptCheck')
    dataSet.TorqueOptCheck = 1;
    dataSet.TorRipOptCheck = 1;
    if Dflag
        disp('rev260 - added check for optimization objective')
    end
    flag = 1;
end

if ~isfield(dataSet,'Qs')
    Q = 6*geo.p*geo.q;
    t2 = gcd(round(dataSet.NumOfSlots*6*dataSet.NumOfPolePairs),2*dataSet.NumOfPolePairs);
    dataSet.Qs = Q/t2;
    clear t2;
    if Dflag
        disp('rev260 - added number of simulated stator slot')
    end
    flag = 1;
end

if ~isfield(dataSet,'xRange')
    dataSet.xRange = [0.5 0.8];
    dataSet.bRange = [0.3 0.7];
    if Dflag
        disp('rev260 - added (x,b) range for syrmDesign')
    end
    flag = 1;
end

if ~isfield(dataSet,'BarFillFac')
    dataSet.BarFillFac = 0;
    if Dflag
        disp('rev260 - added barrier filling factor')
    end
    flag = 1;
end



%% from rev 260 to 261
% disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')
% disp('update to 261')

if ~isfield(dataSet,'MassCuOptCheck')
    dataSet.MassCuOptCheck = 0;
    dataSet.MaxCuMass = 0;
    if Dflag
        disp('rev261 - added copper mass in optimization objective')
    end
    flag = 1;
end

if ~isfield(dataSet,'SlotLayerPosCheck')
    dataSet.SlotLayerPosCheck = 0;
    if Dflag
        disp('rev261 - added slot layer position in GUI')
    end
    flag = 1;
end

if ~isfield(dataSet,'RadRibCheck')
    dataSet.RadRibCheck = 0;
    dataSet.RadRibEdit = zeros(1,geo.nlay);
    if Dflag
        disp('rev261 - added radial ribs in GUI')
    end
    flag = 1;
end

%% Bfe and kt
if ~isfield(dataSet,'Bfe')
    dataSet.Bfe=1.5;
    dataSet.kt=1;
    if Dflag
        disp('rev289 - added Bfe and kt')
    end
    flag = 1;
end

%% PMs area in Seg and ISeg for bonded-to-real PMs (W. Ventura Master Thesis)
if ~isfield(dataSet,'Areavert')
    dataSet.Areavert  = zeros(1,4);
    dataSet.Areaob    = zeros(1,4);
    dataSet.Areatot   = zeros(1,4);
    dataSet.Areavert0 = zeros(1,4);
    dataSet.Areaob0   = zeros(1,4);
    dataSet.dob       = ones(1,4);
    dataSet.dvert     = ones(1,4);
    if Dflag
        disp('rev309 - added PMs area evaluation for Seg and ISeg')
    end
    flag = 1;
end

% In some case, Areavert and its friends are set as empty, but is an error.
if isempty(dataSet.Areavert)
    dataSet.Areavert  = zeros(1,4);
    dataSet.Areaob    = zeros(1,4);
    dataSet.Areatot   = zeros(1,4);
    dataSet.Areavert0 = zeros(1,4);
    dataSet.Areaob0   = zeros(1,4);
    dataSet.dob       = ones(1,4);
    dataSet.dvert     = ones(1,4);
end

%% Multi 3-phase option (Simone Adamo Master Thesis)
if ~isfield(dataSet,'Num3PhaseCircuit') %AS
    dataSet.Num3PhaseCircuit=1;
    dataSet.WinFlag=[1 1 1];
    dataSet.DefaultWinMatr=dataSet.WinMatr;
    if Dflag
        disp('rev319 - added multi 3-phase windings')
    end
    flag = 1;
end

%% Added parameters in dataSet for Vtype rotor geometry - rev.Gallo 17/01/2018 
if ~isfield(dataSet,'SlopeBarrier')

    dataSet.SlopeBarrier=60;
    dataSet.SlopeBarrBou=[10 89];
    dataSet.SlopeBarrBouCheck=0;
    
    if Dflag
        disp('rev325 - added Vtype rotor geometry')
    end
    flag = 1;
end

%% Added parameters in dataSet for optimization Vtype rotor - rev.Gallo 20/03/2018
if ~isfield(dataSet,'MaxPMMass')
    dataSet.MaxPMMass= 1.58; %Massa volume magnete macchina di riferimento (mot.232 Pareto OUT_20180308T163427)
    dataSet.MassPMOptCheck=0;
    if Dflag
        disp('rev325 - added PM Mass in optimization objective')
    end
    flag = 1;
end

%% Split radial ribs
if ~isfield(dataSet,'RadRibSplit')
    dataSet.RadRibSplit=0;
    if Dflag
        disp('rev326 - added splitted radial ribs for Seg geometry')
    end
    flag=1;
end

%% FEAfix syrmDesign
if ~isfield(dataSet,'FEAfixN')
    dataSet.FEAfixN=1;
    if Dflag
        disp('rev326 - added FEAfix for syrmDesign tool')
    end
    flag=1;
end

%% Flux Barrier Shift 
if ~isfield(dataSet,'thetaFBS')
    dataSet.thetaFBS=0;
    dataSet.ThetaFBSBouCheck=0;
    dataSet.ThetaFBSBou=[0 360/(6*dataSet.NumOfPolePairs*dataSet.NumOfSlots)];
    if Dflag
        disp('rev326 - added Flux Barrier Shift for Circular and Seg geometries')
    end
    flag=1;
end

%% BrDesign
if ~isfield(dataSet,'BrDesign')
    matTmp=material_properties_layer('Bonded-Magnet');
    dataSet.BrDesign=matTmp.Br;
    if Dflag
        disp('rev326 - added BrDesign useful for PMs design')
    end
    flag=1;
end

%% rewriting geo, per and mat (and check if mat exist)
% if ~exist('mat')
%     flag_mat=1;
% else
%     flag_mat=0;
% end
[bounds, objs, geo, per, mat] = data0(dataSet);
geo.RQ=buildDefaultRQ(bounds);
% if flag_mat && Dflag
%     disp('rev261 - added mat structure for material properties')
% end

%% check if RQ and RQnames are correct
if (length(dataSet.RQ)~=length(dataSet.RQnames) || length(geo.RQ)~=length(dataSet.RQ))
    dataSet.RQ=buildDefaultRQ(bounds);
    dataSet.RQnames=geo.RQnames;
    if Dflag
        disp('rev274 - correct RQ and RQnames')
    end
    flag=1;
end

%% message in command window if some data are added
if flag && Dflag
    disp('The selected project is in an older version of SyR-e: save machine to update to the latest version')
end

