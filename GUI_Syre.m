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

function varargout = GUI_Syre(varargin)
% GUI_SYRE MATLAB code for GUI_Syre.fig
%      GUI_SYRE,by itself, creates a new GUI_SYRE or raises the existing
%      singleton*.
%
%      H = GUI_SYRE returns the handle to a new GUI_SYRE or the handle to
%      the existing singleton*.
%
%      GUI_SYRE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SYRE.M with the given input arguments.
%
%      GUI_SYRE('Property','Value',...) creates a new GUI_SYRE or raises the
%         existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Syre_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Syre_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Syre

% Last Modified by GUIDE v2.5 08-Mar-2017 13:25:40

% Begin initialization code - DO NOT EDIT
% clc

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_Syre_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_Syre_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);

if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI_Syre is made visible.
function GUI_Syre_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Syre (see VARARGIN)

addpath('./mfiles');

% Set the colors indicating a selected/unselected tab
handles.unselectedTabColor = 0.9*[1  1  1];
handles.selectedTabColor = [1  1  1];

% Choose default command line output for GUI_Syre
handles.output = hObject;

set(handles.GeometricTab,'BackgroundColor',handles.selectedTabColor);
set(handles.StatorTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.MaterialsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.WindingsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.OptionsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.OptimizationTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.PostProcTab,'BackgroundColor',handles.unselectedTabColor);

handles.pan1pos = get(handles.OptimizationPanel,'Position');
handles.left_right_corner = handles.pan1pos(2) + handles.pan1pos(4);  % height
handles.pan2pos = get(handles.GeoPanel,'Position');
handles.pan3pos = get(handles.StatorPanel,'Position');
handles.pan4pos = get(handles.WindingsPanel,'Position');
handles.pan5pos = get(handles.OptionsPanel,'Position');
handles.pan6pos = get(handles.MaterialPanel,'Position');
handles.pan7pos = get(handles.PostProcePanel,'Position');
handles.pan8pos = get(handles.RotorGeometryPanel,'Position');
set(handles.GeoPanel,'Position',[handles.pan1pos(1) (handles.left_right_corner-handles.pan2pos(4)) handles.pan2pos(3) handles.pan2pos(4)])
set(handles.StatorPanel,'Position',[handles.pan1pos(1) (handles.left_right_corner-handles.pan3pos(4)) handles.pan3pos(3) handles.pan3pos(4)])
set(handles.WindingsPanel,'Position',[handles.pan1pos(1) (handles.left_right_corner-handles.pan4pos(4)) handles.pan4pos(3) handles.pan4pos(4)])
set(handles.OptionsPanel,'Position',[handles.pan1pos(1) (handles.left_right_corner-handles.pan5pos(4)) handles.pan5pos(3) handles.pan5pos(4)])
set(handles.MaterialPanel,'Position',[handles.pan1pos(1) (handles.left_right_corner-handles.pan6pos(4)) handles.pan6pos(3) handles.pan6pos(4)])
set(handles.PostProcePanel,'Position',[handles.pan1pos(1) (handles.left_right_corner-handles.pan7pos(4)) handles.pan7pos(3) handles.pan7pos(4)])
set(handles.RotorGeometryPanel,'Position',[(handles.pan1pos(1) + handles.pan3pos(3)) (handles.left_right_corner-handles.pan8pos(4)) handles.pan8pos(3) handles.pan8pos(4)])

set(handles.GeoPanel,'Visible','on')
set(handles.StatorPanel,'Visible','off')
set(handles.WindingsPanel,'Visible','off')
set(handles.OptionsPanel,'Visible','off')
set(handles.MaterialPanel,'Visible','off')
set(handles.OptimizationPanel,'Visible','off')
set(handles.PostProcePanel,'Visible','off')
set(handles.RotorGeometryPanel,'Visible','off')

%% === VISIBLE & ENABLE ===================================================
% set(handles.text22,'Visible','off');
% set(handles.PitchWindEdit,'Visible','off');

set(handles.AlphadegreeEdit,'Enable','off');
set(handles.hcmmEdit,'Enable','off');
set(handles.EstimatedCoppTemp,'Enable','off');
set(handles.CalculatedRatedCurrent,'Enable','off');

%% === FLAG VARI ==========================================================
handles.MatrixWinFlag = 1;
handles.UpdateData0 = 0;
handles.Opti = 0;
%% ========================================================================


% Update handles structure
% guidata(hObject, handles);
axes(handles.axes4);
SyreImg = imread('syre.png');
image(SyreImg);
axis off;
axes(handles.axes5);
set(handles.axes5,'YTickLabel',[]);
set(handles.axes5,'XTickLabel',[]);
box on;

load mot_01.mat

[dataSet,geo,per] = back_compatibility(dataSet,geo,per); % backward compatibility with previous version, to delete here
dataSet.currentpathname = [cd '\'];
dataSet.currentfilename = 'mot_01.mat';

%% ========================================================================
%% update new added parameters of ver. 259 & 260 
% if ~isfield(dataSet,'TorqueOptCheck')
%     dataSet.TorqueOptCheck = 0;
%     dataSet.TorRipOptCheck = 0;
% end
% 
% if ~isfield(dataSet,'LossEvaluationCheck')
%     dataSet.LossEvaluationCheck = 0;
%     dataSet.HysteresisLossFactor = 0;
%     dataSet.HysteresisFrequencyFactor = 0;
%     dataSet.HysteresisFluxDenFactor = 0;
%     dataSet.EddyCurLossFactor = 0;
%     dataSet.IronMassDen = 0;
%     dataSet.EvalSpeed = 0;
% end
handles.dataSet = dataSet; % data structure
SetParameters(handles,dataSet) % aux. function for set the values in the edit boxes
handles = DrawPush_Callback(hObject, eventdata, handles);

%% ======= Matrix of winding visual =======================================
columnName = cell(1,floor(geo.Qs));
for i = 1 : floor(geo.Qs)
    columnName{i} = ['Slot n° ',num2str(i)];
end
rowName{1} = 'Layer 1';
rowName{2} = 'Layer 2';
set(handles.WinTable,'rowname',rowName);
set(handles.WinTable,'columnname',columnName);
set(handles.WinTable,'data',dataSet.WinMatr(:,1:floor(geo.Qs)));
per.tempcuest = temp_est_simpleMod(geo,per);
dataSet.EstimatedCopperTemp = per.tempcuest;
set(handles.EstimatedCoppTemp,'String',num2str(per.tempcuest));
handles.dataSet = dataSet; % data structure
%% ========================================================================
guidata(hObject,handles);

% --- Outputs from this function are returned to the command line.
function varargout = GUI_Syre_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function PolePairsEdit_Callback(hObject, eventdata, handles)
% Number of pole pairs(geo.p)
% hObject    handle to PolePairsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of PolePairsEdit as text
%        str2double(get(hObject,'String')) returns contents of PolePairsEdit as a double
dataSet = handles.dataSet;
dataSet.NumOfPolePairs = str2double(get(hObject,'String'));
p = dataSet.NumOfPolePairs;
Q = round(dataSet.NumOfSlots*6*p);
yq = dataSet.PitchShortFac*dataSet.NumOfSlots*3;
path = pwd;
cd(fullfile (path,'koil'));
system(['koil_syre.exe',' ',num2str(Q),' ',num2str(p),' ',num2str(yq)]);
cd(path);
Windings = MatrixWin();
%dataSet.WinMatr = Windings; % winding matrix
% t = gcd(round(dataSet.NumOfSlots*6*dataSet.NumOfPolePairs),dataSet.NumOfPolePairs);  % periodicity
% if ((6*t/Q)>1)
%     Qs = Q/t;   % periodic machine
% else
%     Qs = Q/2/t; % anti-periodic machine
% end
t2 = gcd(round(dataSet.NumOfSlots*6*dataSet.NumOfPolePairs),2*dataSet.NumOfPolePairs);
Qs = Q/t2;
set(handles.SlotSimulEdit,'String',int2str(Qs));  % set Qs in the GUI with the default value
dataSet.Qs = Qs;
dataSet.WinMatr = Windings(:,1:floor(Qs)); % winding matrix, only Qs columns
columnName = cell(1,floor(Qs));
for i = 1 : floor(Qs)
    columnName{i} = ['Slot n° ',num2str(i)];
end
rowName{1} = 'Layer 1';
rowName{2} = 'Layer 2';
set(handles.WinTable,'rowname',rowName);
set(handles.WinTable,'columnname',columnName);
set(handles.WinTable,'data',dataSet.WinMatr(:,1:floor(Qs)));
dataSet.ToothLength = dataSet.StatorOuterRadius * (1-dataSet.AirGapRadius/dataSet.StatorOuterRadius*(1+dataSet.MagLoadingYoke/dataSet.NumOfPolePairs));
dataSet.ToothLength = roundn(dataSet.ToothLength,-4);
set(handles.ToothLengEdit,'String',mat2str(dataSet.ToothLength));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function PolePairsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PolePairsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% round to a set number of decimals cfr
function num = rrtd(num, cfr)
error(nargchk(2, 2, nargin, 'struct'))
validateattributes(num, {'single', 'double'}, {}, 'ROUNDN', 'NUM')
validateattributes(cfr, ...
    {'numeric'}, {'scalar', 'real', 'integer'}, 'ROUNDN', 'CFR')
if cfr < 0
    cc = 10 ^ -cfr;
    num = round(cc * num) / cc;
elseif cfr > 0
    cc = 10 ^ cfr;
    num = cc * round(num / cc);
else
    num = round(num);
end

%% AUX FUNCTIONS ==========================================================
%% SET PARAMETERS =========================================================
function SetParameters(handles,dataIn)
set(handles.PolePairsEdit,'String',num2str(dataIn.NumOfPolePairs));
set(handles.GapThiEdit,'String',num2str(dataIn.AirGapThickness));
set(handles.StatorOuterRadEdit,'String',num2str(dataIn.StatorOuterRadius));
set(handles.AirGapRadiusEdit,'String',num2str(dataIn.AirGapRadius));
set(handles.ShaftRadEdit,'String',num2str(dataIn.ShaftRadius));
set(handles.StackLenghtEdit,'String',num2str(dataIn.StackLength));
set(handles.NumOfSlotsEdit,'String',num2str(dataIn.NumOfSlots));
set(handles.ToothLengEdit,'String',num2str(dataIn.ToothLength));
set(handles.StatorSlotOpeEdit,'String',num2str(dataIn.StatorSlotOpen));
set(handles.ToothWidthEdit,'String',num2str(dataIn.ToothWidth));
set(handles.ToothTanDepEdit,'String',num2str(dataIn.ToothTangDepth));
set(handles.ToothTangAngleEdit,'String',num2str(dataIn.ToothTangAngle));
set(handles.FillCorSlotEdit,'String',num2str(dataIn.FilletCorner));
set(handles.SlotMatEdit,'String',dataIn.SlotMaterial);
set(handles.StatorMatEdit,'String',dataIn.StatorMaterial);
set(handles.RotorMatEdit,'String',dataIn.RotorMaterial);
set(handles.FluxBarMatEdit,'String',dataIn.FluxBarrierMaterial);
set(handles.ShaftMatEdit,'String',dataIn.ShaftMaterial);
set(handles.RotorConMatEdit,'String',dataIn.RotorCondMaterial);
set(handles.SlotFillFacEdit,'String',num2str(dataIn.SlotFillFactor));
set(handles.PitchWindEdit,'String',num2str(dataIn.PitchShortFac));
set(handles.TurnsSeriesEdit,'String',num2str(dataIn.TurnsInSeries));
set(handles.JouleLossesEdit,'String',num2str(dataIn.AdmiJouleLosses));
set(handles.CopperTempEdit,'String',num2str(dataIn.TargetCopperTemp));
set(handles.HousingTempEdit,'String',num2str(dataIn.HousingTemp));
set(handles.EstimatedCoppTemp,'String',num2str(dataIn.EstimatedCopperTemp));
% set(handles.DCLinkVolEdit,'String',num2str(dataIn.DCVoltage));
set(handles.CurrentOverLoadEdit,'String',num2str(dataIn.CurrOverLoad));
set(handles.NumberOfLayersEdit,'String',num2str(dataIn.NumOfLayers));
set(handles.OverSpeedEdit,'String',num2str(dataIn.OverSpeed));
set(handles.OverSpeedEdit,'String',num2str(dataIn.OverSpeed));
set(handles.BrPMEdit,'String',mat2str(dataIn.Br));
if strcmp(dataIn.FluxBarrierMaterial,'Bonded-Magnet')
    set(handles.BrPMEdit,'Enable','on');
else
    set(handles.BrPMEdit,'Enable','off');
end
set(handles.BrPPEdit,'String',mat2str(dataIn.BrPP));
set(handles.SlotSimulEdit,'String',num2str(dataIn.Qs));

set(handles.BarFillFacEdit,'String',mat2str(dataIn.BarFillFac));

set(handles.MaxGenEdit,'String',num2str(dataIn.MaxGen));
set(handles.XPopEdit,'String',num2str(dataIn.XPop));
set(handles.MinExpTorEdit,'String',num2str(dataIn.MinExpTorque));
set(handles.NGridPPEdit,'String',num2str(dataIn.NumGrid));
set(handles.MaxExpeRippleTorEdit,'String',num2str(dataIn.MaxRippleTorque));
set(handles.MecTolerEdit,'String',num2str(dataIn.MinMechTol));
set(handles.SimPosMOEdit,'String',num2str(dataIn.SimPoMOOA));
set(handles.RotorPosiMOEdit,'String',num2str(dataIn.RotPoMOOA));
set(handles.SimPosFinerEdit,'String',num2str(dataIn.SimPoFine));
set(handles.RotorPosFinerEdit,'String',num2str(dataIn.RotPoFine));
set(handles.MeshEdit,'String',num2str(dataIn.Mesh));
set(handles.MeshMOOAEdit,'String',num2str(dataIn.Mesh_MOOA));
set(handles.NumOfRotorPosiPPEdit,'String',num2str(dataIn.NumOfRotPosPP));
set(handles.SpanEltPPEdit,'String',num2str(dataIn.AngularSpanPP));
set(handles.Alpha1BouEdit,'String',mat2str(dataIn.Alpha1Bou));
set(handles.DeltaAlphaBouEdit,'String',mat2str(dataIn.DeltaAlphaBou));
set(handles.hcBouEdit,'String',mat2str(dataIn.hcBou));
set(handles.DfeBouEdit,'String',mat2str(dataIn.DfeBou));
set(handles.PhaseAngleCurrBouEdit,'String',mat2str(dataIn.PhaseAngleCurrBou));
set(handles.GapBouEdit,'String',mat2str(dataIn.GapBou));
set(handles.BrBouEdit,'String',mat2str(dataIn.BrBou));
set(handles.AirgapRadiusBouEdit,'String',mat2str(dataIn.GapRadiusBou));
set(handles.ToothWidthBouEdit,'String',mat2str(dataIn.ToothWiBou));
set(handles.ToothLenBouEdit,'String',mat2str(dataIn.ToothLeBou));
set(handles.StatorSlotOpenBouEdit,'String',mat2str(dataIn.StatorSlotOpenBou));
set(handles.ToothTangDepthBouEdit,'String',mat2str(dataIn.ToothTangDepthBou));
set(handles.LossEvaluationCheck,'Value',dataIn.LossEvaluationCheck);

if dataIn.LossEvaluationCheck == 0
    set(handles.HysteresisLossFactorEdit,'Enable','off');
    set(handles.HysteresisFrequencyFactorEdit,'Enable','off');
    set(handles.HysteresisFluxDenEdit,'Enable','off');
    set(handles.EddyCurLossFactorEdit,'Enable','off');
    set(handles.MassDensityEdit,'Enable','off');
    set(handles.EvaluatedSpeedEdit,'Enable','off');
    dataIn.HysteresisLossFactor = 0;
    dataIn.HysteresisFrequencyFactor = 0;
    dataIn.HysteresisFluxDenFactor = 0;
    dataIn.EddyCurLossFactor = 0;
    set(handles.HysteresisLossFactorEdit,'String',mat2str(dataIn.HysteresisLossFactor));
    set(handles.HysteresisFrequencyFactorEdit,'String',mat2str(dataIn.HysteresisFrequencyFactor));
    set(handles.HysteresisFluxDenEdit,'String',mat2str(dataIn.HysteresisFluxDenFactor));
    set(handles.EddyCurLossFactorEdit,'String',mat2str(dataIn.EddyCurLossFactor));
    set(handles.MassDensityEdit,'String',mat2str(0));
    set(handles.EvaluatedSpeedEdit,'String',mat2str(0));
else
    set(handles.HysteresisLossFactorEdit,'Enable','on');
    set(handles.HysteresisFrequencyFactorEdit,'Enable','on');
    set(handles.HysteresisFluxDenEdit,'Enable','on');
    set(handles.EddyCurLossFactorEdit,'Enable','on');
    set(handles.MassDensityEdit,'Enable','on');
    set(handles.EvaluatedSpeedEdit,'Enable','on');
    set(handles.HysteresisLossFactorEdit,'String',mat2str(dataIn.HysteresisLossFactor));
    set(handles.HysteresisFrequencyFactorEdit,'String',mat2str(dataIn.HysteresisFrequencyFactor));
    set(handles.HysteresisFluxDenEdit,'String',mat2str(dataIn.HysteresisFluxDenFactor));
    set(handles.EddyCurLossFactorEdit,'String',mat2str(dataIn.EddyCurLossFactor));
    set(handles.MassDensityEdit,'String',mat2str(dataIn.IronMassDen));
    set(handles.EvaluatedSpeedEdit,'String',mat2str(dataIn.EvalSpeed));
%     Q = round(dataIn.NumOfSlots*6*dataIn.NumOfPolePairs);
%     t2 = gcd(Q,2*dataIn.NumOfPolePairs);
%     QsCalc = Q/t2;
%     dataIn.AngularSpanPP = max(60,20*QsCalc/dataIn.NumOfSlots);
%     set(handles.SpanEltPPEdit,'String',mat2str(dataIn.AngularSpanPP));
%     dataIn.NumOfRotPosPP = ceil(dataIn.AngularSpanPP/60)*6;
%     set(handles.NumOfRotorPosiPPEdit,'String',mat2str(dataIn.NumOfRotPosPP));
end

if isfield(dataIn,'bRange')
    set(handles.bRangeEdit,'String',mat2str(dataIn.bRange));
    set(handles.xRangeEdit,'String',mat2str(dataIn.xRange));
else
    set(handles.bRangeEdit,'String','NaN');
    set(handles.xRangeEdit,'String','NaN');
end
set(handles.BfeEdit,'String',num2str(dataIn.Bfe));
set(handles.ktEdit,'String',num2str(dataIn.kt));
set(handles.CurrLoXBEdit,'String',num2str(dataIn.CurrLoPP));
set(handles.CurrLoPPEdit,'String',mat2str(dataIn.CurrLoPP));

set(handles.GammaPPEdit,'String',mat2str(dataIn.GammaPP));
set(handles.DxEdit,'String',mat2str(dataIn.DepthOfBarrier));
% set(handles.RQPlotEdit,'String',mat2str(dataIn.RQ));
set(handles.AlphapuEdit,'String',mat2str(dataIn.ALPHApu));
set(handles.hcpuEdit,'String',mat2str(dataIn.HCpu));

set(handles.currentMotFileName,'String',dataIn.currentfilename);
set(handles.MassCuEdit,'String',dataIn.MaxCuMass);
set(handles.MassCuOptCheck,'Value',dataIn.MassCuOptCheck);
if dataIn.NumOfSlots<1
    set(handles.SlotLayerPosCheck,'Enable','on');
else
    set(handles.SlotLayerPosCheck,'Enable','off');
end
set(handles.SlotLayerPosCheck,'Value',dataIn.SlotLayerPosCheck);
set(handles.RadRibCheck,'Value',dataIn.RadRibCheck);
%set(handles.RadRibEdit,'String',mat2str(dataIn.RadRibEdit));
if dataIn.RadRibCheck == 1
    set(handles.RadRibEdit,'Enable','on');
else
    set(handles.RadRibEdit,'Enable','off');
end

%% CHECK BOU ==============================================================
set(handles.Dalpha1BouCheck,'Value',dataIn.Dalpha1BouCheck);
set(handles.DalphaBouCheck,'Value',dataIn.DalphaBouCheck);
set(handles.hcBouCheck,'Value',dataIn.hcBouCheck);
set(handles.DxBouCheck,'Value',dataIn.DxBouCheck);
set(handles.GammaBouCheck,'Value',dataIn.GammaBouCheck);
set(handles.GapBouCheck,'Value',dataIn.GapBouCheck);
set(handles.BrBouCheck,'Value',dataIn.BrBouCheck);
set(handles.AirgapRadiusBouCheck,'Value',dataIn.AirgapRadiusBouCheck);
set(handles.ToothWidthBouCheck,'Value',dataIn.ToothWidthBouCheck);
set(handles.ToothLengthBouCheck,'Value',dataIn.ToothLengthBouCheck);
set(handles.StatorSlotOpenBouCheck,'Value',dataIn.StatorSlotOpenBouCheck);
set(handles.ToothTangDepthBouCheck,'Value',dataIn.ToothTangDepthBouCheck);
%% ========================================================================

%% CHECK OBJ ==============================================================
set(handles.TorqueOptCheck,'Value',dataIn.TorqueOptCheck);
set(handles.TorRipOptCheck,'Value',dataIn.TorRipOptCheck);
if isfield(dataIn,'MassCuOptCheck')                             % da eliminare con compatibilità
    set(handles.MassCuOptCheck,'Value',dataIn.MassCuOptCheck);
else
    set(handles.MassCuOptCheck,'Value',0);
end
%% ========================================================================

if strcmp(dataIn.RMVTmp,'OFF')
    set(handles.RemTMPRadio,'Value',0);
else
    set(handles.RemTMPRadio,'Value',1);
end

if strcmp(dataIn.XFEMMOpt,'N')
    set(handles.XFEMMOptRadio,'Value',0);
else
    set(handles.XFEMMOptRadio,'Value',1);
end

if strcmp(dataIn.XFEMMPPMot,'N')
    set(handles.PostProcXFEMMCheck,'Value',0);
else
    set(handles.PostProcXFEMMCheck,'Value',1);
end

% Listbox materiali
contents = get(handles.TypeOfRotorList,'String');
N = length(contents);
idx = 1;
for k = 1 : N
    if strcmp(dataIn.TypeOfRotor,contents{k})
        idx = k;
        break;
    end
end
set(handles.TypeOfRotorList,'Value',idx);

if strcmp(dataIn.TypeOfRotor,'SPM')
    [~, ~, ~, ~, mat] = data0(dataIn);
    dataIn.Br = mat.LayerMag.Br;
    set(handles.BrPMEdit,'String',num2str(dataIn.Br));
    set(handles.BrPPEdit,'String',num2str(dataIn.Br));
    dataIn.BrPP = dataIn.Br;
    set(handles.BrPMEdit,'Enable','off');
%% limit b range by 3<lm/g<8
%     dataIn.bRange(1)=3;
%     dataIn.bRange(2)=8;
    dataIn.bRange = round(dataIn.bRange*100)/100;
    set(handles.bRangeEdit,'String',mat2str(dataIn.bRange));
end

handles.dataSet = dataIn;
%% === Values for the edit text Plot ======================================
% [bounds, objs, geo, per, mat] = data0(dataIn);
% if not(isempty(bounds))
%     if sum(strcmp(dataIn.RQnames,'dalpha')) >= 1
%         dataIn.Dalpha1BouCheck = 1;
%         dataIn.DalphaBouCheck = 1;
%     end
%     if sum(strcmp(dataIn.RQnames,'hc')) >= 1
%         dataIn.hcBouCheck = 1;
%     end
%     if sum(strcmp(dataIn.RQnames,'Br')) >= 1
%         dataIn.BrBouCheck = 1;
%     end
%     if sum(strcmp(dataIn.RQnames,'g')) >= 1
%         dataIn.GapBouCheck = 1;
%     end
%     if sum(strcmp(dataIn.RQnames,'r')) >= 1
%         dataIn.AirgapRadiusBouCheck = 1;
%     end
%     if sum(strcmp(dataIn.RQnames,'wt')) >= 1
%         dataIn.ToothWidthBouCheck = 1;
%     end
%     if sum(strcmp(dataIn.RQnames,'lt')) >= 1
%         dataIn.ToothLengthBouCheck = 1;
%     end
%     if sum(strcmp(dataIn.RQnames,'acs')) >= 1
%         dataIn.StatorSlotOpenBouCheck = 1;
%     end
%     if sum(strcmp(dataIn.RQnames,'ttd')) >= 1
%         dataIn.ToothTangDepthBouCheck = 1;
%     end
%     if sum(strcmp(dataIn.RQnames,'gamma')) >= 1
%         dataIn.ToothTangDepthBouCheck = 1;
%     end
% end

%% === SPM ===================================================
if strcmp(dataIn.TypeOfRotor,'SPM')
    set(handles.AlphadegreeEdit,'Enable','on');
    set(handles.hcmmEdit,'Enable','on');
    set(handles.hcpuEdit,'String',mat2str(0));
    set(handles.hcpuEdit,'Enable','off');
    set(handles.AlphapuEdit,'String',mat2str(0));
    set(handles.AlphapuEdit,'Enable','off');
    set(handles.DxEdit,'String',mat2str(dataIn.DepthOfBarrier));
    set(handles.DxEdit,'Enable','on');

    set(handles.NumberOfLayersEdit,'String',mat2str(dataIn.AngleSpanOfPM));
    set(handles.NumberOfLayersEdit,'Enable','off');
    set(handles.AlphadegreeEdit,'String',mat2str(dataIn.AngleSpanOfPM));
    set(handles.hcmmEdit,'String',mat2str(dataIn.ThicknessOfPM));
    set(handles.Alpha1BouEdit,'Enable','off');
    set(handles.Dalpha1BouCheck,'Enable','off','Value',0);
    set(handles.DeltaAlphaBouEdit,'Enable','off');
    set(handles.DalphaBouCheck,'Enable','off','Value',0);
    set(handles.DfeBouEdit,'Enable','off');
    set(handles.DxBouCheck,'Enable','off','Value',0);
    set(handles.RadRibCheck,'Enable','off');
    set(handles.RadRibCheck,'Value',0);
    set(handles.RadRibEdit,'Enable','off');
else
    set(handles.AlphadegreeEdit,'Enable','off');
    set(handles.AlphapuEdit,'Enable','on');
    set(handles.hcmmEdit,'Enable','off');
    set(handles.hcpuEdit,'Enable','on');
    if strcmp(dataIn.TypeOfRotor,'ISeg')
        set(handles.DxEdit,'Enable','off');
        set(handles.DxBouCheck,'Enable','off');
        set(handles.DxBouCheck,'Value',0);
        set(handles.DfeBouEdit,'Enable','off');
    else
        set(handles.DxEdit,'Enable','on');
        set(handles.DxBouCheck,'Enable','on');
        set(handles.DxBouCheck,'Value',dataIn.DxBouCheck);
        set(handles.DfeBouEdit,'Enable','on');
    end
    set(handles.NumberOfLayersEdit,'Enable','on');
    if dataIn.Dalpha1BouCheck
        set(handles.Alpha1BouEdit,'Enable','on');
    end
    set(handles.Alpha1BouEdit,'Enable','on');
    set(handles.Dalpha1BouCheck,'Enable','on','Value',dataIn.Dalpha1BouCheck);
    set(handles.DeltaAlphaBouEdit,'Enable','on');
    set(handles.DalphaBouCheck,'Enable','on','Value',dataIn.DalphaBouCheck);
    
    set(handles.RadRibCheck,'Enable','on');
    set(handles.RadRibCheck,'Value',dataIn.RadRibCheck);
    if dataIn.RadRibCheck
        set(handles.RadRibEdit,'Enable','on');
    else
        set(handles.RadRibEdit,'Enable','off');
    end
end

function GapThiEdit_Callback(hObject, eventdata, handles)
% Airgap Thickness[mm](geo.g)
% hObject    handle to GapThiEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of GapThiEdit as text
%        str2double(get(hObject,'String')) returns contents of GapThiEdit as a double
dataSet = handles.dataSet;
dataSet.AirGapThickness = str2double(get(hObject,'String'));
if strcmp(dataSet.TypeOfRotor,'SPM')
%     dataSet.bRange(1)=3;
%     dataSet.bRange(2)=8;
    dataSet.bRange = round(dataSet.bRange*100)/100;
    set(handles.bRangeEdit,'String',mat2str(dataSet.bRange));
end
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function GapThiEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GapThiEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function StatorOuterRadEdit_Callback(hObject, eventdata, handles)
% Stator outer radius[mm](geo.R)
% hObject    handle to StatorOuterRadEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of StatorOuterRadEdit as text
%        str2double(get(hObject,'String')) returns contents of StatorOuterRadEdit as a double
dataSet = handles.dataSet;
dataSet.StatorOuterRadius = str2double(get(hObject,'String'));
dataSet.ToothLength = dataSet.StatorOuterRadius * (1-dataSet.AirGapRadius/dataSet.StatorOuterRadius*(1+dataSet.MagLoadingYoke/dataSet.NumOfPolePairs));
dataSet.ToothLength = roundn(dataSet.ToothLength,-4);
set(handles.ToothLengEdit,'String',mat2str(dataSet.ToothLength));
set(handles.ToothLenBouEdit,'String',mat2str(dataSet.ToothLeBou));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function StatorOuterRadEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StatorOuterRadEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function AirGapRadiusEdit_Callback(hObject, eventdata, handles)
% Airgap Radius [mm](geo.r)
% hObject    handle to AirGapRadiusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of AirGapRadiusEdit as text
%        str2double(get(hObject,'String')) returns contents of AirGapRadiusEdit as a double
dataSet = handles.dataSet;
dataSet.AirGapRadius = str2double(get(hObject,'String'));
dataSet.ToothLength = dataSet.StatorOuterRadius * (1-dataSet.AirGapRadius/dataSet.StatorOuterRadius*(1+dataSet.MagLoadingYoke/dataSet.NumOfPolePairs));
dataSet.ToothLength = roundn(dataSet.ToothLength,-4);
set(handles.ToothLengEdit,'String',mat2str(dataSet.ToothLength));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function AirGapRadiusEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AirGapRadiusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ShaftRadEdit_Callback(hObject, eventdata, handles)
% Shaft radius [mm] (geo.Ar)
% hObject    handle to ShaftRadEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of ShaftRadEdit as text
%        str2double(get(hObject,'String')) returns contents of ShaftRadEdit as a double
dataSet = handles.dataSet;
dataSet.ShaftRadius = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ShaftRadEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ShaftRadEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function StackLenghtEdit_Callback(hObject, eventdata, handles)
% Stack length [mm] (geo.l)
% hObject    handle to StackLenghtEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of StackLenghtEdit as text
%        str2double(get(hObject,'String')) returns contents of StackLenghtEdit as a double
dataSet = handles.dataSet;
dataSet.StackLength = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function StackLenghtEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StackLenghtEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in TypeOfRotorList.
function TypeOfRotorList_Callback(hObject, eventdata, handles)
% Type of rotor (geo.RotType)
% hObject    handle to TypeOfRotorList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns TypeOfRotorList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TypeOfRotorList
contents = get(hObject,'String');
s = contents{get(hObject,'Value')};
dataSet = handles.dataSet;
dataSet.TypeOfRotor = s;
if (strcmp(dataSet.TypeOfRotor,'Fluid') || strcmp(dataSet.TypeOfRotor,'Seg') || strcmp(dataSet.TypeOfRotor,'Circular'))
    dataSet.DxBouCheck = 0;
    set(handles.DxBouCheck,'Value',dataSet.DxBouCheck);
    set(handles.DxEdit,'Enable','on');
    set(handles.DfeBouEdit,'Enable','on');
else
    dataSet.DxBouCheck = 0;
    set(handles.DxBouCheck,'Value',dataSet.DxBouCheck);
    set(handles.DxEdit,'Enable','off');
    set(handles.DfeBouEdit,'Enable','off');
end
if strcmp(dataSet.TypeOfRotor,'SPM')
    %set(handles.AlphadegreeEdit,'Enable','on');
    %set(handles.hcmmEdit,'Enable','on');
    %set(handles.hcpuEdit,'String',mat2str(0));
    %set(handles.hcpuEdit,'Enable','off');
    %set(handles.AlphapuEdit,'String',mat2str(0));
    %set(handles.AlphapuEdit,'Enable','off');
    %set(handles.NumberOfLayersEdit,'String',mat2str(0));
    %set(handles.NumberOfLayersEdit,'Enable','off');
    %set(handles.AlphadegreeEdit,'String',mat2str(dataSet.AngleSpanOfPM));
    %set(handles.hcmmEdit,'String',mat2str(dataSet.ThicknessOfPM));
    %set(handles.DxEdit,'Enable','on');
    dataSet.DepthOfBarrier = 1;
    %set(handles.DxEdit,'String',mat2str(dataSet.DepthOfBarrier)); 
    %set(handles.Alpha1BouEdit,'Enable','off');
    %set(handles.DeltaAlphaBouEdit,'Enable','off');
    dataSet.Dalpha1BouCheck =0;
    dataSet.DalphaBouCheck = 0;
    dataSet.DxBouCheck = 0;
    %set(handles.DfeBouEdit,'Enable','off');
    dataSet.FluxBarrierMaterial = 'BMN-38EH';
    [~,~,~,~,mat] = data0(dataSet);
    dataSet.Br = mat.LayerMag.Br;
    dataSet.BrPP = dataSet.Br;
    set(handles.FluxBarMatEdit,'String',dataSet.FluxBarrierMaterial);
    set(handles.BrPMEdit,'String',num2str(dataSet.Br));
    set(handles.BrPPEdit,'String',num2str(dataSet.BrPP));
    %set(handles.RadRibCheck,'Value',0);
    %set(handles.RadRibEdit,'Enable','off');
    %   3<lm/g<8
%     dataSet.bRange(1)=3;
%     dataSet.bRange(2)=8;
    dataSet.bRange = [4 6];
    dataSet.bRange = round(dataSet.bRange*100)/100;
    set(handles.bRangeEdit,'String',mat2str(dataSet.bRange));
%     [bounds, objs, geo, per] = data0(dataSet);
else
    %set(handles.AlphadegreeEdit,'Enable','off');        
    %set(handles.hcmmEdit,'Enable','off');
    %set(handles.hcpuEdit,'Enable','on');
    set(handles.hcpuEdit,'String',mat2str(dataSet.HCpu));
    set(handles.AlphapuEdit,'Enable','on');
    set(handles.AlphapuEdit,'String',mat2str(dataSet.ALPHApu));
    set(handles.NumberOfLayersEdit,'Enable','on');
    set(handles.NumberOfLayersEdit,'String',mat2str(dataSet.NumOfLayers));
    set(handles.Alpha1BouEdit,'Enable','on');
    set(handles.DeltaAlphaBouEdit,'Enable','on');
    if (strcmp(dataSet.TypeOfRotor,'Fluid') || strcmp(dataSet.TypeOfRotor,'Seg') || strcmp(dataSet.TypeOfRotor,'Circular'))
        if not(length(dataSet.DepthOfBarrier) == dataSet.NumOfLayers)
            dataSet.DepthOfBarrier = zeros(1,dataSet.NumOfLayers);
        end
        set(handles.DxEdit,'String',mat2str(dataSet.DepthOfBarrier));
    end
    
    dataSet.FluxBarrierMaterial = 'Air';
    [bounds, objs, geo, per, mat] = data0(dataSet);
    data = buildDefaultRQ(bounds);
    dataSet.RQ = data;
    dataSet.Br = mat.LayerMag.Br;
    dataSet.BrPP = dataSet.Br;
    set(handles.FluxBarMatEdit,'String',dataSet.FluxBarrierMaterial);
    set(handles.BrPMEdit,'String',num2str(dataSet.Br));
    set(handles.BrPPEdit,'String',num2str(dataSet.BrPP));
    if max(dataSet.bRange)>1
        dataSet.bRange = [0.4 0.6];
    end
    
    
    
%     if (strcmp(dataSet.TypeOfRotor,'Fluid') || strcmp(dataSet.TypeOfRotor,'Seg'))
%         set(handles.DxEdit,'String',mat2str(dataSet.DepthOfBarrier));
%     end
    [truefalse, index] = ismember('dx', geo.RQnames);
    if truefalse
        dataSet.DfeBou = bounds(index,:);            % barrier offset [p.u.]
    end
    flag_plot = 'Y';
    h = handles.axes5;
    [hc,dalpha] = Plot_Machine(h,dataSet,flag_plot);
    view = round(100*[dalpha hc])/100;
    set(handles.AlphadegreeEdit,'String',mat2str(view(1:dataSet.NumOfLayers)));
    set(handles.hcmmEdit,'String',mat2str(view((dataSet.NumOfLayers+1):end)));
    set(handles.DfeBouEdit,'String',mat2str(dataSet.DfeBou));
    set(handles.bRangeEdit,'String',mat2str(dataSet.bRange));
end
handles.dataSet = dataSet;
SetParameters(handles,dataSet);
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function TypeOfRotorList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TypeOfRotorList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function NumOfSlotsEdit_Callback(hObject, eventdata, handles)
% Number of slots/(poles*phases) (geo.q)
% hObject    handle to NumOfSlotsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of NumOfSlotsEdit as text
%        str2double(get(hObject,'String')) returns contents of NumOfSlotsEdit as a double
dataSet = handles.dataSet;
dataSet.NumOfSlots = eval(get(hObject,'String'));
% dataSet.yq = dataSet.NumOfSlots*3;
% set(handles.YqEdit,'String',num2str(dataSet.yq));
p = dataSet.NumOfPolePairs;
Q = round(dataSet.NumOfSlots*6*p);
yq = dataSet.PitchShortFac*dataSet.NumOfSlots*3;
path = pwd;
cd(fullfile (path,'koil'));
system(['koil_syre.exe',' ',num2str(Q),' ',num2str(p),' ',num2str(yq)]);
cd(path);
Windings = MatrixWin();
%dataSet.WinMatr = Windings; % winding matrix
% t = gcd(round(dataSet.NumOfSlots*6*dataSet.NumOfPolePairs),dataSet.NumOfPolePairs);  % periodicity
% if ((6*t/Q)>1)
%     Qs = Q/t;   % periodic machine
% else
%     Qs = Q/2/t; % anti-periodic machine
% end
t2 = gcd(round(dataSet.NumOfSlots*6*dataSet.NumOfPolePairs),2*dataSet.NumOfPolePairs);
Qs = Q/t2;
set(handles.SlotSimulEdit,'String',int2str(Qs));  % set Qs in the GUI with the default value
dataSet.Qs = Qs;
dataSet.WinMatr = Windings(:,1:floor(Qs)); % winding matrix, only Qs columns
columnName = cell(1,floor(Qs));
for i = 1 : floor(Qs)
    columnName{i} = ['Slot n° ',num2str(i)];
end
rowName{1} = 'Layer 1';
rowName{2} = 'Layer 2';
set(handles.WinTable,'rowname',rowName);
set(handles.WinTable,'columnname',columnName);
set(handles.WinTable,'data',dataSet.WinMatr(:,1:floor(Qs)));
if dataSet.NumOfSlots<1
    set(handles.SlotLayerPosCheck,'Enable','on');
else
    set(handles.SlotLayerPosCheck,'Enable','off');
end
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function NumOfSlotsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumOfSlotsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ToothLengEdit_Callback(hObject, eventdata, handles)
% Tooth length [mm](geo.lt)
% hObject    handle to ToothLengEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of ToothLengEdit as text
%        str2double(get(hObject,'String')) returns contents of ToothLengEdit as a double
dataSet = handles.dataSet;
dataSet.ToothLength = str2double(get(hObject,'String'));
if dataSet.ToothLength>(dataSet.StatorOuterRadius - dataSet.AirGapRadius-dataSet.AirGapThickness)
    dataSet.ToothLength = dataSet.StatorOuterRadius - dataSet.AirGapRadius-dataSet.AirGapThickness;
    set(handles.ToothLengEdit,'String',mat2str(dataSet.ToothLength));
end
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ToothLengEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ToothLengEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function StatorSlotOpeEdit_Callback(hObject, eventdata, handles)
% Stator slot opening [p.u.] (geo.acs)
% hObject    handle to StatorSlotOpeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of StatorSlotOpeEdit as text
%        str2double(get(hObject,'String')) returns contents of StatorSlotOpeEdit as a double
dataSet = handles.dataSet;
dataSet.StatorSlotOpen = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function StatorSlotOpeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StatorSlotOpeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ToothWidthEdit_Callback(hObject, eventdata, handles)
% Tooth width [mm](geo.wt)
% hObject    handle to ToothWidthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of ToothWidthEdit as text
%        str2double(get(hObject,'String')) returns contents of ToothWidthEdit as a double
dataSet = handles.dataSet;
dataSet.ToothWidth = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ToothWidthEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ToothWidthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ToothTanDepEdit_Callback(hObject, eventdata, handles)
% Tooth tang. depth [mm](geo.ttd)
% hObject    handle to ToothTanDepEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of ToothTanDepEdit as text
%        str2double(get(hObject,'String')) returns contents of ToothTanDepEdit as a double
dataSet = handles.dataSet;
dataSet.ToothTangDepth = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ToothTanDepEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ToothTanDepEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ToothTangAngleEdit_Callback(hObject, eventdata, handles)
% Tooth tang. angle [mech. deg.](geo.tta)
% hObject    handle to ToothTangAngleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of ToothTangAngleEdit as text
%        str2double(get(hObject,'String')) returns contents of ToothTangAngleEdit as a double
dataSet = handles.dataSet;
dataSet.ToothTangAngle = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ToothTangAngleEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ToothTangAngleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FillCorSlotEdit_Callback(hObject, eventdata, handles)
% Fillet at the back iron corner of the slot[mm](geo.SFR)
% hObject    handle to FillCorSlotEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of FillCorSlotEdit as text
%        str2double(get(hObject,'String')) returns contents of FillCorSlotEdit as a double
dataSet = handles.dataSet;
dataSet.FilletCorner = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function FillCorSlotEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FillCorSlotEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SlotFillFacEdit_Callback(hObject, eventdata, handles)
% Slot filling factor(geo.Kcu)
% hObject    handle to SlotFillFacEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of SlotFillFacEdit as text
%        str2double(get(hObject,'String')) returns contents of SlotFillFacEdit as a double
dataSet = handles.dataSet;
dataSet.SlotFillFactor = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function SlotFillFacEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SlotFillFacEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PitchWindEdit_Callback(hObject, eventdata, handles)
% Pitch shortening factor(geo.Kracc)
% hObject    handle to PitchWindEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of PitchWindEdit as text
%        str2double(get(hObject,'String')) returns contents of PitchWindEdit as a double
dataSet = handles.dataSet;
dataSet.PitchShortFac = eval(get(hObject,'String'));
p = dataSet.NumOfPolePairs;
Q = round(dataSet.NumOfSlots*6*p);
yq = dataSet.PitchShortFac*dataSet.NumOfSlots*3;
path = pwd;
cd(fullfile (path,'koil'));
system(['koil_syre.exe',' ',num2str(Q),' ',num2str(p),' ',num2str(yq)]);
cd(path);
Windings = MatrixWin();
%dataSet.WinMatr = Windings; % winding matrix
% t = gcd(round(dataSet.NumOfSlots*6*dataSet.NumOfPolePairs),dataSet.NumOfPolePairs);  % periodicity
% if ((6*t/Q)>1)
%     Qs = Q/t;   % periodic machine
% else
%     Qs = Q/2/t; % anti-periodic machine
% end
t2 = gcd(round(dataSet.NumOfSlots*6*dataSet.NumOfPolePairs),2*dataSet.NumOfPolePairs);
Qs = Q/t2;
set(handles.SlotSimulEdit,'String',int2str(Qs));  % set Qs in the GUI with the default value
dataSet.Qs = Qs;
dataSet.WinMatr = Windings(:,1:floor(Qs)); % winding matrix, only Qs columns
columnName = cell(1,floor(Qs));
for i = 1 : floor(Qs)
    columnName{i} = ['Slot n° ',num2str(i)];
end
rowName{1} = 'Layer 1';
rowName{2} = 'Layer 2';
set(handles.WinTable,'rowname',rowName);
set(handles.WinTable,'columnname',columnName);
set(handles.WinTable,'data',dataSet.WinMatr(:,1:floor(Qs)));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function PitchWindEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PitchWindEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TurnsSeriesEdit_Callback(hObject, eventdata, handles)
% Turns in series per phase(geo.Ns)
% hObject    handle to TurnsSeriesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of TurnsSeriesEdit as text
%        str2double(get(hObject,'String')) returns contents of TurnsSeriesEdit as a double
dataSet = handles.dataSet;
dataSet.TurnsInSeries = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function TurnsSeriesEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TurnsSeriesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function JouleLossesEdit_Callback(hObject, eventdata, handles)
% permitted Joule losses at stall (per.Loss)
% hObject    handle to JouleLossesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of JouleLossesEdit as text
%        str2double(get(hObject,'String')) returns contents of JouleLossesEdit as a double
dataSet = handles.dataSet;
dataSet.AdmiJouleLosses = str2double(get(hObject,'String'));
% [bounds, objs, geo, per, mat] = data0(dataSet);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function JouleLossesEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to JouleLossesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CopperTempEdit_Callback(hObject, eventdata, handles)
% copper temperature
% hObject    handle to CopperTempEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of CopperTempEdit as text
%        str2double(get(hObject,'String')) returns contents of CopperTempEdit as a double
dataSet = handles.dataSet;
dataSet.TargetCopperTemp = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function CopperTempEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CopperTempEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% function DCLinkVolEdit_Callback(hObject, eventdata, handles)
% % hObject    handle to DCLinkVolEdit (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % Hints: get(hObject,'String') returns contents of DCLinkVolEdit as text
% %        str2double(get(hObject,'String')) returns contents of DCLinkVolEdit as a double
% dataSet = handles.dataSet;
% dataSet.DCVoltage = str2double(get(hObject,'String'));
% handles.dataSet = dataSet;
% guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function DCLinkVolEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DCLinkVolEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CurrentOverLoadEdit_Callback(hObject, eventdata, handles)
% current overload factor used for optimization (per overload)
% hObject    handle to CurrentOverLoadEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of CurrentOverLoadEdit as text
%        str2double(get(hObject,'String')) returns contents of CurrentOverLoadEdit as a double
dataSet = handles.dataSet;
dataSet.CurrOverLoad = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function CurrentOverLoadEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrentOverLoadEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function NumberOfLayersEdit_Callback(hObject, eventdata, handles)
% Number of rotor barriers(geo.nlay)
% hObject    handle to NumberOfLayersEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of NumberOfLayersEdit as text
%        str2double(get(hObject,'String')) returns contents of NumberOfLayersEdit as a double
dataSet = handles.dataSet;
dataSet.NumOfLayers = str2double(get(hObject,'String'));
% if (dataSet.NumOfLayers==1)
%     dataSet.DalphaBouCheck = 0;
%     set(handles.DalphaBouCheck,'Value',dataSet.DalphaBouCheck);
%     set(handles.DeltaAlphaBouEdit,'Enable','off');
% else
%     dataSet.DalphaBouCheck = 1;
%     set(handles.DalphaBouCheck,'Value',dataSet.DalphaBouCheck);
%     set(handles.DeltaAlphaBouEdit,'Enable','on');
% end
[bounds, objs, geo, per, mat] = data0(dataSet);
data = buildDefaultRQ(bounds);
dataSet.RQ = data;
if(dataSet.DalphaBouCheck)
    last_index = dataSet.NumOfLayers;
    dataSet.ALPHApu = data(1:last_index);
else
    last_index = 0;
    dataSet.ALPHApu = ones(1,dataSet.NumOfLayers)*roundn(1/(dataSet.NumOfLayers+0.5),-2);
end
if(dataSet.hcBouCheck)
    dataSet.HCpu = data((last_index+1):end-1);
else
    dataSet.HCpu = ones(1,dataSet.NumOfLayers)*0.5;
end
dataSet.HCpu = round(dataSet.HCpu,2);
dataSet.DepthOfBarrier = zeros(1,dataSet.NumOfLayers);
set(handles.AlphapuEdit,'String',mat2str(dataSet.ALPHApu));
set(handles.hcpuEdit,'String',mat2str(dataSet.HCpu));
set(handles.DxEdit,'String',mat2str(dataSet.DepthOfBarrier));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function NumberOfLayersEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumberOfLayersEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function OverSpeedEdit_Callback(hObject, eventdata, handles)
% overspeed (geo.nmax)
% hObject    handle to OverSpeedEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of OverSpeedEdit as text
%        str2double(get(hObject,'String')) returns contents of OverSpeedEdit as a double
dataSet = handles.dataSet;
dataSet.OverSpeed = str2double(get(hObject,'String'));
dataSet.RadRibCheck = 0;
dataSet.RadRibEdit = zeros(1,dataSet.NumOfLayers);
handles.dataSet = dataSet;
SetParameters(handles,dataSet)
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function OverSpeedEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OverSpeedEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function BrPMEdit_Callback(hObject, eventdata, handles)
% Br for PMs
% hObject    handle to BrPMEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of BrPMEdit as text
%        str2double(get(hObject,'String')) returns contents of BrPMEdit as a double
dataSet = handles.dataSet;
% dataSet.Br = str2double(get(hObject,'String'));
dataSet.Br = str2num(get(hObject,'String'));
dataSet.BrPP = dataSet.Br;  % update Br post proc when Br is changed

set(handles.BrPMEdit,'String',mat2str(dataSet.Br));
set(handles.BrPPEdit,'String',mat2str(dataSet.BrPP));

dataSet.BrBouCheck = 0;
set(handles.BrBouCheck,'Value',dataSet.BrBouCheck);

% SetParameters(handles,dataSet);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function BrPMEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BrPMEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function HcPMEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HcPMEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function ScalingFacEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScalingFacEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function MecTolerEdit_Callback(hObject, eventdata, handles)
% minimum mechanical tolerance (geo.pont0)
% hObject    handle to MecTolerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of MecTolerEdit as text
%        str2double(get(hObject,'String')) returns contents of MecTolerEdit as a double
dataSet = handles.dataSet;
dataSet.MinMechTol = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function MecTolerEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MecTolerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SimPosFinerEdit_Callback(hObject, eventdata, handles)
% Rotor position for Pareto(geo.nsim_singt)
% hObject    handle to SimPosFinerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of SimPosFinerEdit as text
%        str2double(get(hObject,'String')) returns contents of SimPosFinerEdit as a double
dataSet = handles.dataSet;
dataSet.SimPoFine = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function SimPosFinerEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SimPosFinerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function RotorPosFinerEdit_Callback(hObject, eventdata, handles)
% Rotor angular excursion for Pareto
% hObject    handle to RotorPosFinerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of RotorPosFinerEdit as text
%        str2double(get(hObject,'String')) returns contents of RotorPosFinerEdit as a double
dataSet = handles.dataSet;
dataSet.RotPoFine = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function RotorPosFinerEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RotorPosFinerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SimPosMOEdit_Callback(hObject, eventdata, handles)
% Simuluated position for MODE(geo.nsim_MOOA)
% hObject    handle to SimPosMOEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of SimPosMOEdit as text
%        str2double(get(hObject,'String')) returns contents of SimPosMOEdit as a double
dataSet = handles.dataSet;
dataSet.SimPoMOOA = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function SimPosMOEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SimPosMOEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function RotorPosiMOEdit_Callback(hObject, eventdata, handles)
% Rotor angular excursion for MODE (geo.delta_sim_MOOA)
% hObject    handle to RotorPosiMOEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of RotorPosiMOEdit as text
%        str2double(get(hObject,'String')) returns contents of RotorPosiMOEdit as a double
dataSet = handles.dataSet;
dataSet.RotPoMOOA = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function RotorPosiMOEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RotorPosiMOEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SlotMatEdit_Callback(hObject, eventdata, handles)
% Slot material(geo.BLKLABELSmaterials)
% hObject    handle to SlotMatEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of SlotMatEdit as text
%        str2double(get(hObject,'String')) returns contents of SlotMatEdit as a double
dataSet = handles.dataSet;
dataSet.SlotMaterial = get(hObject,'String');
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function SlotMatEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SlotMatEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function StatorMatEdit_Callback(hObject, eventdata, handles)
% stator material(geo.BLKLABELSmaterials)
% hObject    handle to StatorMatEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of StatorMatEdit as text
%        str2double(get(hObject,'String')) returns contents of StatorMatEdit as a double
dataSet = handles.dataSet;
dataSet.StatorMaterial = get(hObject,'String');
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function StatorMatEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StatorMatEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function RotorMatEdit_Callback(hObject, eventdata, handles)
% rotor material (geo.BLKLABELSmaterials)
% hObject    handle to RotorMatEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of RotorMatEdit as text
%        str2double(get(hObject,'String')) returns contents of RotorMatEdit as a double
dataSet = handles.dataSet;
dataSet.RotorMaterial = get(hObject,'String');
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function RotorMatEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RotorMatEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FluxBarMatEdit_Callback(hObject, eventdata, handles)
% flux barrier material (geo.BLKLABELSmaterials)
% hObject    handle to FluxBarMatEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of FluxBarMatEdit as text
%        str2double(get(hObject,'String')) returns contents of FluxBarMatEdit as a double
dataSet = handles.dataSet;
dataSet.FluxBarrierMaterial = get(hObject,'String');
mat = material_properties_layer(dataSet.FluxBarrierMaterial);
dataSet.Br = mat.Br;
dataSet.BrPP = dataSet.Br;
if strcmp(dataSet.TypeOfRotor,'SPM')
%% limit b range by 3<lm/g<8
%     dataSet.bRange(1)=3;
%     dataSet.bRange(2)=8;
    dataSet.bRange = round(dataSet.bRange*100)/100;
    set(handles.bRangeEdit,'String',mat2str(dataSet.bRange));
end

handles.dataSet = dataSet;
SetParameters(handles,dataSet);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function FluxBarMatEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FluxBarMatEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ShaftMatEdit_Callback(hObject, eventdata, handles)
% shaft material (geo.BLKLABELSmaterials)
% hObject    handle to ShaftMatEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of ShaftMatEdit as text
%        str2double(get(hObject,'String')) returns contents of ShaftMatEdit as a double
dataSet = handles.dataSet;
dataSet.ShaftMaterial = get(hObject,'String');
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ShaftMatEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ShaftMatEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function RotorConMatEdit_Callback(hObject, eventdata, handles)
% rotor conductor winding material (geo.BLKLABELSmaterials)
% hObject    handle to RotorConMatEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of RotorConMatEdit as text
%        str2double(get(hObject,'String')) returns contents of RotorConMatEdit as a double
dataSet = handles.dataSet;
dataSet.RotorCondMaterial = get(hObject,'String');
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function RotorConMatEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RotorConMatEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in LibraryMaterPush.
% function LibraryMaterPush_Callback(hObject, eventdata, handles)
% % Open Library Material
% % hObject    handle to LibraryMaterPush (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% fileID = fopen('empty_case.fem');
% C = textscan(fileID,'%s');
% fclose(fileID);
% n = length(C{1});
% A = zeros(100,2);
% j = 0;
% for i = 3 : 1 : n
%     if strcmp(C{1}(i-2),'<BlockName>') && strcmp(C{1}(i-1),'=')
%         j = j + 1;
%         A(j,1) = i; %matrice degli indici dei materiali
%         while ~strcmp(C{1}(i+1),'<Mu_x>')
%             i = i + 1;
%         end
%         A(j,2) = i;
%     end
% end
% [m,~] = size(A);
% for i = 1 : m
%     if A(i,:) == 0
%         m = i - 1;
%         break
%     end
% end
% material = {};
% for i = 1 : m
%     if A(i,1) == A(i,2)
%         material = strjoin([material,C{1}(A(i))]);
%     else
%         diff = A(i,2) - A(i,1);
%         index = A(i,1);
%         name = {};
%         for j = 0 : diff
%             name = strcat([name,C{1}(index + j)]);
%         end
%         material = strjoin([material,name]);
%     end
% end
% 
% set(handles.MaterialText,'Style','Edit');
% set(handles.MaterialText,'Max',m);
% set(handles.MaterialText,'String',material);
% guidata(hObject,handles)

% --- Executes on button press in OptimizePush.
function OptimizePush_Callback(hObject, eventdata, handles)
% Optimize pushbutton
% hObject    handle to OptimizePush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataSet = handles.dataSet;
save('dataSet','dataSet');
figure()
[bounds, objs, geo, per, mat] = data0(dataSet);
dat.geo0 = geo;
per.objs = objs;
dat.per = per;
dat.mat = mat;
%%%%%%%%%% FEMM fitness handle %%%%%%%%%%%%%%%%%%%%%%%%%%
eval_type = 'MO_OA';
if strcmp(dataSet.XFEMMOpt,'Y')
    FitnessFunction = @(x)FEMMfitnessX(x,geo,per,eval_type);
else
    FitnessFunction = @(x)FEMMfitness(x,geo,per,mat,eval_type);
end
dat.CostProblem = FitnessFunction;           % Cost function instance
%% ========================================================================
NOBJ = length(objs);
NOBJ = sum(objs(:,2));
XPOP = dataSet.XPop;
Esc = 0.75;
Pm= 0.2;
NVAR = size(bounds,1);
MAXGEN = dataSet.MaxGen;
MAXFUNEVALS = 20000*NVAR*NOBJ;
%% Variables regarding the optimization problem
dat.FieldD = bounds;
dat.Initial = bounds;
dat.NOBJ = NOBJ;
dat.NRES = 0;
dat.NVAR = NVAR;
dat.mop = str2func('evaluateF');
dat.eval_type = eval_type;
dat.XPOP = XPOP;
dat.Esc = Esc;
dat.Pm  = Pm;
dat.fl  = 0.1;
dat.fu  = 0.9;
dat.tau1= 0.1;
dat.tau2= 0.1;
dat.InitialPop=[];
dat.MAXGEN = MAXGEN;
dat.MAXFUNEVALS = MAXFUNEVALS;
dat.SaveResults='yes';
dat.CounterGEN=0;
dat.CounterFES=0;

% Run the algorithm.
OUT = MODE2(dat, dataSet);

delete('dataSet.mat');
rmdir('partial_optimization','s');

% guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function YqEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YqEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in UpdateDataPush.
function UpdateDataPush_Callback(hObject, eventdata, handles)
% hObject    handle to UpdateDataPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
flag = 1;
% save([cd,'\tmp\GUI\flag.mat'],'flag');
save([cd,'\tmp\flag.mat'],'flag');
dataSet = handles.dataSet;
% save([cd,'\tmp\GUI\dataSet.mat'],'dataSet');
save([cd,'\tmp\dataSet.mat'],'dataSet');
guidata(hObject,handles)


% --------------------------------------------------------------------
% function SavePushTool_ClickedCallback(hObject, eventdata, handles)
% % hObject    handle to SavePushTool (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% dataSet = handles.dataSet;
% uisave('dataSet','Project.mat');
% guidata(hObject,handles);

% --------------------------------------------------------------------
function LoadPushTool_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to LoadPushTool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.mat');
load([PathName FileName]);
if ~exist('geo')
    geo = geo0;
end

if exist('dataSet')~=1
    dataSet=build_dataSet(geo,per);
    disp('dataSet reconstructed. Please check the data')
end
[dataSet,geo,per] = back_compatibility(dataSet,geo,per);
dataSet.RQ = roundn(dataSet.RQ,-4);
GeometricTab_Callback(hObject, eventdata, handles)

dataSet.currentpathname = PathName;
dataSet.currentfilename = FileName;

%% updata new added data of version 259 & 260
% if ~isfield(dataSet,'TorqueOptCheck')
%     dataSet.TorqueOptCheck = 0;
%     dataSet.TorRipOptCheck = 0;
% end
% 
% if ~isfield(dataSet,'LossEvaluationCheck')
%     dataSet.LossEvaluationCheck = 0;
%     dataSet.HysteresisLossFactor = 0;
%     dataSet.HysteresisFrequencyFactor = 0;
%     dataSet.HysteresisFluxDenFactor = 0;
%     dataSet.EddyCurLossFactor = 0;
%     dataSet.IronMassDen = 0;
%     dataSet.EvalSpeed = 0;
% end
handles.dataSet = dataSet;
SetParameters(handles,dataSet);
handles = DrawPush_Callback(hObject, eventdata, handles);
% handles.dataSet = dataSet;

%% ======= Matrix of winding ==============================================
[~,n] = size(dataSet.WinMatr);
columnName = cell(1,n);
for i = 1 : n
    columnName{i} = ['Slot n° ',num2str(i)];
end
rowName{1} = 'Layer 1';
rowName{2} = 'Layer 2';
set(handles.WinTable,'rowname',rowName);
set(handles.WinTable,'columnname',columnName);
set(handles.WinTable,'data',dataSet.WinMatr);

%% ========================================================================
guidata(hObject,handles);

% --- Executes on button press in GeometricTab.
function GeometricTab_Callback(hObject, eventdata, handles)
% % hObject    handle to GeometricTab (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%
set(handles.GeometricTab,'BackgroundColor',handles.selectedTabColor);
set(handles.StatorTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.MaterialsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.WindingsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.OptionsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.OptimizationTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.PostProcTab,'BackgroundColor',handles.unselectedTabColor);

set(handles.GeoPanel,'Visible','on')
set(handles.StatorPanel,'Visible','off')
set(handles.WindingsPanel,'Visible','off')
set(handles.OptionsPanel,'Visible','off')
set(handles.MaterialPanel,'Visible','off')
set(handles.PostProcePanel,'Visible','off')
set(handles.OptimizationPanel,'Visible','off')
set(handles.RotorGeometryPanel,'Visible','off');


% set(handles.WindingsTab,'Visible','on')
% set(handles.MaterialsTab,'Visible','on')
% set(handles.OptimizationTab,'Visible','on')

guidata(hObject,handles)

% --- Executes on button press in StatorTab.
function StatorTab_Callback(hObject, eventdata, handles)
% hObject    handle to StatorTab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.GeometricTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.StatorTab,'BackgroundColor',handles.selectedTabColor);
set(handles.MaterialsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.WindingsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.OptionsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.OptimizationTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.PostProcTab,'BackgroundColor',handles.unselectedTabColor);

set(handles.GeoPanel,'Visible','off')
set(handles.StatorPanel,'Visible','on')
set(handles.RotorGeometryPanel,'Visible','on');
set(handles.WindingsPanel,'Visible','off')
set(handles.OptionsPanel,'Visible','off')
set(handles.PostProcePanel,'Visible','off')
set(handles.MaterialPanel,'Visible','off')
set(handles.OptimizationPanel,'Visible','off')

% set(handles.WindingsTab,'Visible','on')
% set(handles.MaterialsTab,'Visible','on')
% set(handles.OptimizationTab,'Visible','on')

guidata(hObject,handles)

% --- Executes on button press in MaterialsTab.
function MaterialsTab_Callback(hObject, eventdata, handles)
% hObject    handle to MaterialsTab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.GeometricTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.StatorTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.MaterialsTab,'BackgroundColor',handles.selectedTabColor);
set(handles.WindingsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.OptionsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.OptimizationTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.PostProcTab,'BackgroundColor',handles.unselectedTabColor);

set(handles.GeoPanel,'Visible','off')
set(handles.StatorPanel,'Visible','off')
set(handles.WindingsPanel,'Visible','off')
set(handles.OptionsPanel,'Visible','off')
set(handles.PostProcePanel,'Visible','off')
set(handles.MaterialPanel,'Visible','on')
set(handles.OptimizationPanel,'Visible','off')
set(handles.RotorGeometryPanel,'Visible','off')

guidata(hObject,handles)

% --- Executes on button press in WindingsTab.
function WindingsTab_Callback(hObject, eventdata, handles)
% hObject    handle to WindingsTab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.GeometricTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.StatorTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.MaterialsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.WindingsTab,'BackgroundColor',handles.selectedTabColor);
set(handles.OptionsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.OptimizationTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.PostProcTab,'BackgroundColor',handles.unselectedTabColor);

set(handles.GeoPanel,'Visible','off')
set(handles.StatorPanel,'Visible','off')
set(handles.WindingsPanel,'Visible','on')
set(handles.OptionsPanel,'Visible','off')
set(handles.MaterialPanel,'Visible','off')
set(handles.PostProcePanel,'Visible','off')
set(handles.OptimizationPanel,'Visible','off')
set(handles.RotorGeometryPanel,'Visible','off')

guidata(hObject,handles)

% --- Executes on button press in OptionsTab.
function OptionsTab_Callback(hObject, eventdata, handles)
% hObject    handle to OptionsTab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.GeometricTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.StatorTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.MaterialsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.WindingsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.OptionsTab,'BackgroundColor',handles.selectedTabColor);
set(handles.OptimizationTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.PostProcTab,'BackgroundColor',handles.unselectedTabColor);

set(handles.RotorGeometryPanel,'Visible','off');
set(handles.GeoPanel,'Visible','off')
set(handles.StatorPanel,'Visible','off')
set(handles.WindingsPanel,'Visible','off')
set(handles.OptionsPanel,'Visible','on')
set(handles.MaterialPanel,'Visible','off')
set(handles.PostProcePanel,'Visible','off')
set(handles.OptimizationPanel,'Visible','off')

% set(handles.WindingsTab,'Visible','on')
% set(handles.MaterialsTab,'Visible','on')
% set(handles.OptimizationTab,'Visible','on')

guidata(hObject,handles)

% --- Executes on button press in OptimizationTab.
function OptimizationTab_Callback(hObject, eventdata, handles)
% hObject    handle to OptimizationTab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.GeometricTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.StatorTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.MaterialsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.WindingsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.OptionsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.OptimizationTab,'BackgroundColor',handles.selectedTabColor);
set(handles.PostProcTab,'BackgroundColor',handles.unselectedTabColor);

set(handles.RotorGeometryPanel,'Visible','off');
set(handles.GeoPanel,'Visible','off');
set(handles.StatorPanel,'Visible','off');
set(handles.WindingsPanel,'Visible','off');
set(handles.OptionsPanel,'Visible','off');
set(handles.MaterialPanel,'Visible','off');
set(handles.PostProcePanel,'Visible','off');
set(handles.OptimizationPanel,'Visible','on');

guidata(hObject,handles)

% --- Executes on button press in SaveConfPush.
function SaveConfPush_Callback(hObject, eventdata, handles)
% Save Configuration pushbutton
% hObject    handle to SaveConfPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.WinTable,'data');
dataSet = handles.dataSet;
[m,n] = size(data);
for i = 1 : m
    for k = 1 : n
        if abs(data(i,k)) > 3 || isnan(data(i,k)) || data(i,k) == 0
            dataSet = handles.dataSet;
            if exist([cd,'\tmp\flag.mat']) > 0
                flag = 0;
                save([cd,'\tmp\flag.mat'],'flag');
                flag = 1;
                save([cd,'\tmp\flag.mat'],'flag');
            else
                [~, ~, ~, geo, ~, ~] = data0();
            end
            data = geo.avv;
            msgbox('You must insert valid values','!!Warning!!');
            set(handles.WinTable,'data',data);
            return
        end
    end
end
if m == 1
    data(2,:) = data(1,:);
end
dataSet.WinMatr = data;
msgbox('Matrix of windings saved');
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes on key press with focus on WinTable and none of its controls.
function WinTable_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to WinTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

% --- Executes when entered data in editable cell(s) in WinTable.
function WinTable_CellEditCallback(hObject, eventdata, handles)
% winding layout
% hObject    handle to WinTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in DrawPush.
function handles = DrawPush_Callback(hObject, eventdata, handles)
% hObject    handle to DrawPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
flag_plot = 'Y';
h = handles.axes5;
dataSet = handles.dataSet;
% if ~isfield(dataSet,'TorqueOptCheck')
%     dataSet.TorqueOptCheck = 0;
%     dataSet.TorRipOptCheck = 0;
%     set(handles.TorqueOptCheck,'Value',dataSet.TorqueOptCheck);
%     set(handles.TorRipOptCheck,'Value',dataSet.TorRipOptCheck);
%     handles.dataSet = dataSet;
% end
[~, ~, ~, per, mat] = data0(dataSet);
[hc,dalpha,geo] = Plot_Machine(h,dataSet,flag_plot);

% tempcu, io, alpha, hc refreshed every time display is refreshed
per.tempcuest = temp_est_simpleMod(geo,per);
dataSet.EstimatedCopperTemp = per.tempcuest;
set(handles.EstimatedCoppTemp,'String',num2str(per.tempcuest));
per.io = calc_io(geo,per);
set(handles.CalculatedRatedCurrent,'String',num2str(per.io));
temp = round(100*[dalpha hc])/100;
if ~strcmp (geo.RotType, 'SPM')
    set(handles.AlphadegreeEdit,'String',mat2str(temp(1:dataSet.NumOfLayers)));
    set(handles.hcmmEdit,'String',mat2str(temp((dataSet.NumOfLayers+1):end)));
end
tmp = round(geo.pont*100)/100;
set(handles.RadRibEdit,'String',mat2str(tmp));
dataSet.RadRibEdit = tmp;
dataSet.FilletCorner = geo.SFR;
set(handles.FillCorSlotEdit,'String',mat2str(round(100*geo.SFR)/100));
handles.dataSet = dataSet;
guidata(hObject,handles)

function MeshEdit_Callback(hObject, eventdata, handles)
% mesh (geo.K_mesh)
% hObject    handle to MeshEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of MeshEdit as text
%        str2double(get(hObject,'String')) returns contents of MeshEdit as a double
dataSet = handles.dataSet;
dataSet.Mesh = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function MeshEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MeshEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MeshMOOAEdit_Callback(hObject, eventdata, handles)
% Mesh MOOA (geo.K_mesh_MOOA)
% hObject    handle to MeshMOOAEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of MeshMOOAEdit as text
%        str2double(get(hObject,'String')) returns contents of MeshMOOAEdit as a double
dataSet = handles.dataSet;
dataSet.Mesh_MOOA = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function MeshMOOAEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MeshMOOAEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NumOfRotorPosiPPEdit_Callback(hObject, eventdata, handles)
% Number of rotor positions
% hObject    handle to NumOfRotorPosiPPEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of NumOfRotorPosiPPEdit as text
%        str2double(get(hObject,'String')) returns contents of NumOfRotorPosiPPEdit as a double
dataSet = handles.dataSet;
dataSet.NumOfRotPosPP = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function NumOfRotorPosiPPEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumOfRotorPosiPPEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SpanEltPPEdit_Callback(hObject, eventdata, handles)
% Rotor angular excursion[elect. deg.]
% hObject    handle to SpanEltPPEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of SpanEltPPEdit as text
%        str2double(get(hObject,'String')) returns contents of SpanEltPPEdit as a double
dataSet = handles.dataSet;
dataSet.AngularSpanPP = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function SpanEltPPEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SpanEltPPEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PostProcTab.
function PostProcTab_Callback(hObject, eventdata, handles)
% hObject    handle to PostProcTab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.GeometricTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.StatorTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.MaterialsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.WindingsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.OptionsTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.OptimizationTab,'BackgroundColor',handles.unselectedTabColor);
set(handles.PostProcTab,'BackgroundColor',handles.selectedTabColor);

set(handles.GeoPanel,'Visible','off')
set(handles.StatorPanel,'Visible','off')
set(handles.WindingsPanel,'Visible','off')
set(handles.OptionsPanel,'Visible','off')
set(handles.MaterialPanel,'Visible','off')
set(handles.PostProcePanel,'Visible','on')
set(handles.RotorGeometryPanel,'Visible','off');
set(handles.OptimizationPanel,'Visible','off')

guidata(hObject,handles)



function Alpha1BouEdit_Callback(hObject, eventdata, handles)
% first barrier airgap angle or dalpha 1 (bounds_dalpha_1)
% hObject    handle to Alpha1BouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of Alpha1BouEdit as text
%        str2double(get(hObject,'String')) returns contents of Alpha1BouEdit as a double
dataSet = handles.dataSet;
dataSet.Alpha1Bou = str2num(get(hObject,'String'));
[m,n] = size(dataSet.Alpha1Bou);
if m*n > 2 || m*n == 1
    disp('Must be vector 2x1')
    return
end
[bounds, objs, geo, per, ~] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function Alpha1BouEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Alpha1BouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DeltaAlphaBouEdit_Callback(hObject, eventdata, handles)
% Other barriers angles [p.u.](bounds_dalpha)
% hObject    handle to DeltaAlphaBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of DeltaAlphaBouEdit as text
%        str2double(get(hObject,'String')) returns contents of DeltaAlphaBouEdit as a double
dataSet = handles.dataSet;
dataSet.DeltaAlphaBou = str2num(get(hObject,'String'));
[m,n] = size(dataSet.DeltaAlphaBou);
if m*n > 2 || m*n == 1
    disp('Must be vector 2x1')
    return
end
[bounds, objs, geo, per, ~] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function DeltaAlphaBouEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DeltaAlphaBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hcBouEdit_Callback(hObject, eventdata, handles)
% hc[p.u.](bounds_hc)
% hObject    handle to hcBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of hcBouEdit as text
%        str2double(get(hObject,'String')) returns contents of hcBouEdit as a double
dataSet = handles.dataSet;
dataSet.hcBou = str2num(get(hObject,'String'));
[m,n] = size(dataSet.hcBou);
if m*n > 2 || m*n == 1
    disp('Must be vector 2x1')
    return
end
[bounds, objs, geo, per, ~] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function hcBouEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hcBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DfeBouEdit_Callback(hObject, eventdata, handles)
% dx [p.u.](bounds_Dx)
% hObject    handle to DfeBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of DfeBouEdit as text
%        str2double(get(hObject,'String')) returns contents of DfeBouEdit as a double
dataSet = handles.dataSet;
dataSet.DfeBou = str2num(get(hObject,'String'));
[m,n] = size(dataSet.DfeBou);
if m*n > 2 || m*n == 1
    disp('Must be vector 2x1')
    return
end
[bounds, objs, geo, per, ~] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function DfeBouEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DfeBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function PhaseAngleCurrBouEdit_Callback(hObject, eventdata, handles)
% Current phase angle  [elect. deg.](bounds_gamma)
% hObject    handle to PhaseAngleCurrBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of PhaseAngleCurrBouEdit as text
%        str2double(get(hObject,'String')) returns contents of PhaseAngleCurrBouEdit as a double
dataSet = handles.dataSet;
dataSet.PhaseAngleCurrBou = str2num(get(hObject,'String'));
[m,n] = size(dataSet.PhaseAngleCurrBou);
if m*n > 2 || m*n == 1
    disp('Must be vector 2x1')
    return
end
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function PhaseAngleCurrBouEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PhaseAngleCurrBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% function DeltaXBouEdit_Callback(hObject, eventdata, handles)
% % hObject    handle to DeltaXBouEdit (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % Hints: get(hObject,'String') returns contents of DeltaXBouEdit as text
% %        str2double(get(hObject,'String')) returns contents of DeltaXBouEdit as a double
% dataSet = handles.dataSet;
% dataSet.DeltaXBou = str2num(get(hObject,'String'));
% [m,n] = size(dataSet.DeltaXBou);
% if m*n > 2 || m*n == 1
%     disp('Must be vector 2x1')
%     return
% end
% handles.dataSet = dataSet;
% guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function DeltaXBouEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DeltaXBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RemTMPRadio.
function RemTMPRadio_Callback(hObject, eventdata, handles)
% remove TMP file (geo.RemoveTMPfile)it only works with XFEMM
% hObject    handle to RemTMPRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of RemTMPRadio
dataSet = handles.dataSet;
a = get(hObject,'Value');
if a > 1
    dataSet.RMVTmp = 'ON';
else
    dataSet.RMVTmp = 'OFF';
end
handles.dataSet = dataSet;
guidata(hObject,handles)


function BrPPEdit_Callback(hObject, eventdata, handles)
% Br [T]
% hObject    handle to BrPPEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of BrPPEdit as text
%        str2double(get(hObject,'String')) returns contents of BrPPEdit as a double

dataSet = handles.dataSet;
dataSet.BrPP = str2num((get(hObject,'String')));
set(handles.BrPPEdit,'String',mat2str(dataSet.BrPP));
handles.dataSet = dataSet;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function BrPPEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BrPPEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function GammaPPEdit_Callback(hObject, eventdata, handles)
% Current phase angle[elect. deg.]
% hObject    handle to GammaPPEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of GammaPPEdit as text
%        str2double(get(hObject,'String')) returns contents of GammaPPEdit as a double

dataSet = handles.dataSet;
try
    dataSet.GammaPP = str2num((get(hObject,'String')));
catch
    dataSet.GammaPP = str2double((get(hObject,'String')));
end
handles.dataSet = dataSet;
guidata(hObject,handles)




% --- Executes during object creation, after setting all properties.
function GammaPPEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GammaPPEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CurrLoPPEdit_Callback(hObject, eventdata, handles)
% Current load [p.u.]
% hObject    handle to CurrLoPPEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of CurrLoPPEdit as text
%        str2double(get(hObject,'String')) returns contents of CurrLoPPEdit as a double

dataSet = handles.dataSet;
try
    dataSet.CurrLoPP = str2num((get(hObject,'String')));
catch
    dataSet.CurrLoPP = str2double((get(hObject,'String')));
end
handles.dataSet = dataSet;
set(handles.CurrLoXBEdit,'String',num2str(dataSet.CurrLoPP));
set(handles.CurrLoPPEdit,'String',mat2str(dataSet.CurrLoPP));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function CurrLoPPEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrLoPPEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in StartPProPush.
function StartPProPush_Callback(hObject, eventdata, handles)
% Start pushbutton of processing
% hObject    handle to StartPProPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataSet = handles.dataSet;
post_proc_single_motor(dataSet);
guidata(hObject,handles)


% --- Executes on button press in PostProcXFEMMCheck.
function PostProcXFEMMCheck_Callback(hObject, eventdata, handles)
% Post Processing with XFEMM check box
% hObject    handle to PostProcXFEMMCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of PostProcXFEMMCheck

dataSet = handles.dataSet;
if get(hObject,'Value') > 0
    dataSet.XFEMMPPMot = 'Y';
else
    dataSet.XFEMMPPMot = 'N';
end
handles.dataSet = dataSet;
guidata(hObject,handles)


% --- Executes on button press in XFEMMOptRadio.
function XFEMMOptRadio_Callback(hObject, eventdata, handles)
% XFEMM option check box
% hObject    handle to XFEMMOptRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of XFEMMOptRadio

dataSet = handles.dataSet;
if get(hObject,'Value') > 0
    dataSet.XFEMMOpt = 'Y';
else
    dataSet.XFEMMOpt = 'N';
end
handles.dataSet = dataSet;
guidata(hObject,handles)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over MaterialText.
function MaterialText_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to MaterialText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% function RQPlotEdit_Callback(hObject, eventdata, handles)
% % hObject    handle to RQPlotEdit (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % Hints: get(hObject,'String') returns contents of RQPlotEdit as text
% %        str2double(get(hObject,'String')) returns contents of RQPlotEdit as a double
% dataSet = handles.dataSet;
% dataSet.RQ = str2num(get(hObject,'String'));
% flag_plot = 'Y';
% h = handles.axes5;
% [hc,dalpha] = Plot_Machine(h,dataSet,flag_plot);
% view = round(100*[dalpha hc])/100;
% set(handles.ViewHcDaEdit,'String',mat2str(view));
% handles.dataSet = dataSet;
% guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function RQPlotEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RQPlotEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MatrixOfData0Push.
function MatrixOfData0Push_Callback(hObject, eventdata, handles)
% Matrix of data0
% hObject    handle to MatrixOfData0Push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of MatrixOfData0Push
dataSet = handles.dataSet;
if exist([cd,'\tmp\flag.mat']) > 0
    flag = 0;
    save([cd,'\tmp\flag.mat'],'flag');
    [bounds, ~, geo, ~, ~] = data0();
    flag = 1;
    save([cd,'\tmp\flag.mat'],'flag');
else
    [bounds, ~, geo, ~, ~] = data0();
end
data = geo.avv;
% data = dataSet.WinMatr;
[~,n] = size(data);
columnName = cell(1,n);
for i = 1 : n
    columnName{i} = ['Slot n° ',num2str(i)];
end
rowName{1} = 'Layer 1';
rowName{2} = 'Layer 2';
set(handles.WinTable,'rowname',rowName);
set(handles.WinTable,'columnname',columnName);
set(handles.WinTable,'data',data);
handles.dataSet = dataSet;
guidata(hObject,handles)


function MaxGenEdit_Callback(hObject, eventdata, handles)
% Maximum number of generations
% hObject    handle to MaxGenEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of MaxGenEdit as text
%        str2double(get(hObject,'String')) returns contents of MaxGenEdit as a double
dataSet = handles.dataSet;
dataSet.MaxGen = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function MaxGenEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxGenEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function XPopEdit_Callback(hObject, eventdata, handles)
% Population Size
% hObject    handle to XPopEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of XPopEdit as text
%        str2double(get(hObject,'String')) returns contents of XPopEdit as a double
dataSet = handles.dataSet;
dataSet.XPop = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function XPopEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XPopEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ReturnData0Push.
function ReturnData0Push_Callback(hObject, eventdata, handles)
% hObject    handle to ReturnData0Push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
flag = 0;
save([cd,'\tmp\flag.mat'],'flag');
run GUI_Syre.m
guidata(hObject,handles)


function NGridPPEdit_Callback(hObject, eventdata, handles)
% Number of points in [0 Imax] (n_grid)
% hObject    handle to NGridPPEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of NGridPPEdit as text
%        str2double(get(hObject,'String')) returns contents of NGridPPEdit as a double
dataSet = handles.dataSet;
try
    dataSet.NumGrid = str2num((get(hObject,'String')));
catch
    dataSet.NumGrid = str2double((get(hObject,'String')));
end
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function NGridPPEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NGridPPEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in ClearCachePush.
function ClearCachePush_Callback(hObject, eventdata, handles)
% Empty \tmp folder pushbutton
% hObject    handle to ClearCachePush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    rmdir([cd,'\tmp'],'s');
    msgbox('directory empty');
catch
    msgbox('directory empty');
end
if exist([cd,'\tmp'],'dir') == 0
    mkdir([cd,'\tmp']);
end
guidata(hObject,handles)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% % --- Executes when OptimizationPanel is resized.
% function OptimizationPanel_ResizeFcn(hObject, eventdata, handles)
% % hObject    handle to OptimizationPanel (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)



function ViewHcDaEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ViewHcDaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ViewHcDaEdit as text
%        str2double(get(hObject,'String')) returns contents of ViewHcDaEdit as a double


% --- Executes during object creation, after setting all properties.
function ViewHcDaEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ViewHcDaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SaveMachinePush.
function SaveMachinePush_Callback(hObject, eventdata, handles)
% Save machine pushbutton
% hObject    handle to SaveMachinePush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dataSet = DrawPushMachine(handles);
handles.dataSet = dataSet;
guidata(hObject,handles)


% function DeltaXBouEdit_Callback(hObject, eventdata, handles)
% % hObject    handle to DeltaXBouEdit (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%
% % Hints: get(hObject,'String') returns contents of DeltaXBouEdit as text
% %        str2double(get(hObject,'String')) returns contents of DeltaXBouEdit as a double


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over DrawPush.
% function DrawPush_ButtonDownFcn(hObject, eventdata, handles)
% % hObject    handle to DrawPush (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)



function AlphapuEdit_Callback(hObject, eventdata, handles)
% alpha [p.u.]
% hObject    handle to AlphapuEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hO bject,'String') returns contents of AlphapuEdit as text
%        str2double(get(hObject,'String')) returns contents of AlphapuEdit as a double
dataSet = handles.dataSet;
ALPHApu = str2num(get(hObject,'String'));
% the sum of pu angles is rescaled to one
% alpha1 is not rescaled: only the angles from the second barrier onwards
if sum(ALPHApu)>(1-0.05) % 0.05 is the angular space guaranteed for the spider
    if ALPHApu(1)>(1-0.05*(dataSet.NumOfLayers)) % max ALPHApu(1)=1-0.05-0.05*(nlay-1)
        disp('Wrong value of first alpha_pu')
        ALPHApu(1)=1-0.05*(dataSet.NumOfLayers);
    end
    ALPHApu(2:end)=ALPHApu(2:end)/sum(ALPHApu(2:end))*(1-0.05-ALPHApu(1));
    disp('Correct value of alpha_pu')
end
dataSet.ALPHApu=round(ALPHApu,3);
set(hObject,'String',mat2str(dataSet.ALPHApu));
dataSet.Dalpha1BouCheck = 0;
dataSet.DalphaBouCheck = 0;
set(handles.Dalpha1BouCheck,'Value',dataSet.Dalpha1BouCheck);
set(handles.DalphaBouCheck,'Value',dataSet.DalphaBouCheck);
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function AlphapuEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AlphapuEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function AlphadegreeEdit_Callback(hObject, eventdata, handles)
% alpha [degrees]
% hObject    handle to AlphadegreeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of AlphadegreeEdit as text
%        str2double(get(hObject,'String')) returns contents of AlphadegreeEdit as a double
dataSet = handles.dataSet;
dataSet.AngleSpanOfPM = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
flag_plot = 'Y';
h = handles.axes5;
[hc,dalpha] = Plot_Machine(h,dataSet,flag_plot);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function AlphadegreeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AlphadegreeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function hcpuEdit_Callback(hObject, eventdata, handles)
% hc [p.u.]
% hObject    handle to hcpuEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of hcpuEdit as text
%        str2double(get(hObject,'String')) returns contents of hcpuEdit as a double
dataSet = handles.dataSet;
dataSet.HCpu = str2num(get(hObject,'String'));
dataSet.hcBouCheck = 0;
set(handles.hcBouCheck,'Value',dataSet.hcBouCheck);
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function hcpuEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hcpuEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function hcmmEdit_Callback(hObject, eventdata, handles)
% hc [mm]
% hObject    handle to hcmmEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of hcmmEdit as text
%        str2double(get(hObject,'String')) returns contents of hcmmEdit as a double
dataSet = handles.dataSet;
dataSet.ThicknessOfPM = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
flag_plot = 'Y';
h = handles.axes5;
[hc,dalpha] = Plot_Machine(h,dataSet,flag_plot);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function hcmmEdit_CreateFcn(hObject, ~, handles)
% hObject    handle to hcmmEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BrBouEdit_Callback(hObject, eventdata, handles)
% Br [T](bounds_Br)
% hObject    handle to BrBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BrBouEdit as text
%        str2double(get(hObject,'String')) returns contents of BrBouEdit as a double
dataSet = handles.dataSet;
dataSet.BrBou = str2num(get(hObject,'String'));
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function BrBouEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BrBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BrBouCheck.
function BrBouCheck_Callback(hObject, eventdata, handles)
% Br [T] (bounds_Br)
% hObject    handle to BrBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of BrBouCheck
dataSet = handles.dataSet;
dataSet.BrBouCheck = get(hObject,'Value');
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);

handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)




function AirgapRadiusBouEdit_Callback(hObject, eventdata, handles)
% Airgap radius [mm](bounds_xr)
% hObject    handle to AirgapRadiusBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AirgapRadiusBouEdit as text
%        str2double(get(hObject,'String')) returns contents of AirgapRadiusBouEdit as a double
dataSet = handles.dataSet;
dataSet.GapRadiusBou = str2num(get(hObject,'String'));
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function AirgapRadiusBouEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AirgapRadiusBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AirgapRadiusBouCheck.
function AirgapRadiusBouCheck_Callback(hObject, eventdata, handles)
% Airgap radius [mm](bounds_xr)check box
% hObject    handle to AirgapRadiusBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AirgapRadiusBouCheck
dataSet = handles.dataSet;
dataSet.AirgapRadiusBouCheck = get(hObject,'Value');
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

function ToothWidthBouEdit_Callback(hObject, eventdata, handles)
% Tooth width [mm](bounds_wt)
% hObject    handle to ToothWidthBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ToothWidthBouEdit as text
%        str2double(get(hObject,'String')) returns contents of ToothWidthBouEdit as a double
dataSet = handles.dataSet;
dataSet.ToothWiBou = str2num(get(hObject,'String'));
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ToothWidthBouEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ToothWidthBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ToothWidthBouCheck.
function ToothWidthBouCheck_Callback(hObject, eventdata, handles)
% Tooth width [mm](bounds_wt)check box
% hObject    handle to ToothWidthBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ToothWidthBouCheck
dataSet = handles.dataSet;
dataSet.ToothWidthBouCheck = get(hObject,'Value');
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)


function ToothLenBouEdit_Callback(hObject, eventdata, handles)
% Tooth length [mm](bounds_lt)
% hObject    handle to ToothLenBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ToothLenBouEdit as text
%        str2double(get(hObject,'String')) returns contents of ToothLenBouEdit as a double
dataSet = handles.dataSet;
dataSet.ToothLeBou = str2num(get(hObject,'String'));
if dataSet.ToothLeBou(2)>(dataSet.StatorOuterRadius - dataSet.AirGapRadius-dataSet.AirGapThickness)
    dataSet.ToothLeBou(2) = dataSet.StatorOuterRadius - dataSet.AirGapRadius-dataSet.AirGapThickness;
    set(handles.ToothLenBouEdit,'String',mat2str(dataSet.ToothLeBou));
end
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ToothLenBouEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ToothLenBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ToothLengthBouCheck.
function ToothLengthBouCheck_Callback(hObject, eventdata, handles)
% Tooth length [mm] (bounds_lt) check box
% hObject    handle to ToothLengthBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ToothLengthBouCheck
dataSet = handles.dataSet;
dataSet.ToothLengthBouCheck = get(hObject,'Value');
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)


function GapBouEdit_Callback(hObject, eventdata, handles)
% Airgap thickness[mm](bounds_g)
% hObject    handle to GapBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GapBouEdit as text
%        str2double(get(hObject,'String')) returns contents of GapBouEdit as a double
dataSet = handles.dataSet;
dataSet.GapBou = str2num(get(hObject,'String'));
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function GapBouEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GapBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Dalpha1BouCheck.
function Dalpha1BouCheck_Callback(hObject, eventdata, handles)
% 1st barrier airgap angle [p.u.](bounds_dalpha_1)check box
% hObject    handle to Dalpha1BouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Dalpha1BouCheck
dataSet = handles.dataSet;
dataSet.Dalpha1BouCheck = get(hObject,'Value');
if dataSet.Dalpha1BouCheck ==0
    dataSet.DalphaBouCheck = 0;
else
    dataSet.DalphaBouCheck = 1;
end
set(handles.DalphaBouCheck,'Value',dataSet.DalphaBouCheck);
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes on button press in DalphaBouCheck.
function DalphaBouCheck_Callback(hObject, eventdata, handles)
% Other barriers angles [p.u.](bounds_dalpha)
% hObject    handle to DalphaBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DalphaBouCheck
dataSet = handles.dataSet;
dataSet.DalphaBouCheck = get(hObject,'Value');
if dataSet.DalphaBouCheck ==0
    dataSet.Dalpha1BouCheck = 0;
else
    dataSet.Dalpha1BouCheck = 1;
end
set(handles.Dalpha1BouCheck,'Value',dataSet.Dalpha1BouCheck);
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes on button press in hcBouCheck.
function hcBouCheck_Callback(hObject, eventdata, handles)
% hc [p.u.](bounds_hc) check box
% hObject    handle to hcBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hcBouCheck
dataSet = handles.dataSet;
dataSet.hcBouCheck = get(hObject,'Value');
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes on button press in DxBouCheck.
function DxBouCheck_Callback(hObject, eventdata, handles)
% dx [p.u.](bounds_Dx)
% hObject    handle to DxBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DxBouCheck
dataSet = handles.dataSet;
dataSet.DxBouCheck = get(hObject,'Value');
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)


% --- Executes on button press in GammaBouCheck.
function GammaBouCheck_Callback(hObject, eventdata, handles)
% Current phase angle[elect. deg.](bounds_gamma)
% hObject    handle to GammaBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GammaBouCheck
dataSet = handles.dataSet;
dataSet.GammaBouCheck = get(hObject,'Value');
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes on button press in GapBouCheck.
function GapBouCheck_Callback(hObject, eventdata, handles)
% Airgap thickness[mm](bounds_g)
% hObject    handle to GapBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of GapBouCheck
dataSet = handles.dataSet;
dataSet.GapBouCheck = get(hObject,'Value');
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)


function ToothTangDepthBouEdit_Callback(hObject, eventdata, handles)
% tooth tang depth[mm]  (bounds_ttd)
% hObject    handle to ToothTangDepthBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ToothTangDepthBouEdit as text
%        str2double(get(hObject,'String')) returns contents of ToothTangDepthBouEdit as a double
dataSet = handles.dataSet;
dataSet.ToothTangDepthBou = str2num(get(hObject,'String'));
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function ToothTangDepthBouEdit_CreateFcn(hObject, eventdata, handles)
% tooth tang depth[mm]  (bounds_ttd)
% hObject    handle to ToothTangDepthBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in StatorSlotOpenBouCheck.
function StatorSlotOpenBouCheck_Callback(hObject, eventdata, handles)
% stator slot open [p.u.] (bounds_acs)
% hObject    handle to StatorSlotOpenBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of StatorSlotOpenBouCheck
dataSet = handles.dataSet;
dataSet.StatorSlotOpenBouCheck = get(hObject,'Value');
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)


% --- Executes on button press in ToothTangDepthBouCheck.
function ToothTangDepthBouCheck_Callback(hObject, eventdata, handles)
% tooth tang depth[mm]  (bounds_ttd)
% hObject    handle to ToothTangDepthBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ToothTangDepthBouCheck
dataSet = handles.dataSet;
dataSet.ToothTangDepthBouCheck = get(hObject,'Value');
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)



function StatorSlotOpenBouEdit_Callback(hObject, eventdata, handles)
% stator slot open [p.u.] (bounds_acs)
% hObject    handle to StatorSlotOpenBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StatorSlotOpenBouEdit as text
%        str2double(get(hObject,'String')) returns contents of StatorSlotOpenBouEdit as a double
dataSet = handles.dataSet;
dataSet.StatorSlotOpenBou = str2num(get(hObject,'String'));
[bounds, objs, geo, per, mat] = data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function StatorSlotOpenBouEdit_CreateFcn(hObject, eventdata, handles)
% stator slot open [p.u.] (bounds_acs)
% hObject    handle to StatorSlotOpenBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function EstimatedCoppTemp_Callback(hObject, eventdata, handles)
% hObject    handle to EstimatedCoppTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EstimatedCoppTemp as text
%        str2double(get(hObject,'String')) returns contents of EstimatedCoppTemp as a double


% --- Executes during object creation, after setting all properties.
function EstimatedCoppTemp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EstimatedCoppTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function HousingTempEdit_Callback(hObject, eventdata, handles)
% hObject    handle to HousingTempEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HousingTempEdit as text
%        str2double(get(hObject,'String')) returns contents of HousingTempEdit as a double
dataSet = handles.dataSet;
dataSet.HousingTemp = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function HousingTempEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HousingTempEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DxEdit_Callback(hObject, eventdata, handles)
% hObject    handle to DxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DxEdit as text
%        str2double(get(hObject,'String')) returns contents of DxEdit as a double
dataSet = handles.dataSet;
dataSet.DepthOfBarrier = str2num(get(hObject,'String'));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function DxEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CalculatedRatedCurrent_Callback(hObject, eventdata, handles)
% hObject    handle to CalculatedRatedCurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CalculatedRatedCurrent as text
%        str2double(get(hObject,'String')) returns contents of CalculatedRatedCurrent as a double


% --- Executes during object creation, after setting all properties.
function CalculatedRatedCurrent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CalculatedRatedCurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function currentMotFileName_Callback(hObject, eventdata, handles)
% hObject    handle to currentMotFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currentMotFileName as text
%        str2double(get(hObject,'String')) returns contents of currentMotFileName as a double


% --- Executes during object creation, after setting all properties.
function currentMotFileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentMotFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in syrmDesignPushButt.
function syrmDesignPushButt_Callback(hObject, eventdata, handles)
% hObject    handle to syrmDesignPushButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dataSet = handles.dataSet;
[dataSet,flagS] = syrmDesign(dataSet);
handles.dataSet = dataSet;
if flagS
    dataSet = DrawPushMachine(handles,dataSet.currentfilename,dataSet.currentpathname);
end
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
%flag_plot = 'Y';
%h = handles.axes5;
SetParameters(handles,dataSet)
guidata(hObject,handles)
% [~,~,~] = data0(dataSet);
%[~,~,~] = Plot_Machine(h,dataSet,flag_plot);

function xRangeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to xRangeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xRangeEdit as text
%        str2double(get(hObject,'String')) returns contents of xRangeEdit as a double
dataSet = handles.dataSet;
dataSet.xRange = str2num(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function xRangeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xRangeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bRangeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to bRangeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bRangeEdit as text
%        str2double(get(hObject,'String')) returns contents of bRangeEdit as a double
dataSet = handles.dataSet;
dataSet.bRange = str2num(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function bRangeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bRangeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CurrLoXBEdit_Callback(hObject, eventdata, handles)
% hObject    handle to CurrLoXBEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurrLoXBEdit as text
%        str2double(get(hObject,'String')) returns contents of CurrLoXBEdit as a double
dataSet = handles.dataSet;
try
    dataSet.CurrLoPP = str2num((get(hObject,'String')));
catch
    dataSet.CurrLoPP = str2double((get(hObject,'String')));
end
handles.dataSet = dataSet;
set(handles.CurrLoXBEdit,'String',num2str(dataSet.CurrLoPP));
set(handles.CurrLoPPEdit,'String',mat2str(dataSet.CurrLoPP));
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function CurrLoXBEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrLoXBEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BarFillFacEdit_Callback(hObject, eventdata, handles)
% hObject    handle to BarFillFacEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BarFillFacEdit as text
%        str2double(get(hObject,'String')) returns contents of BarFillFacEdit as a double
dataSet = handles.dataSet;
% dataSet.BarFillFac = str2double(get(hObject,'String'));
dataSet.BarFillFac = str2num(get(hObject,'String'));
if ~strcmp(dataSet.TypeOfRotor,'SPM')
    if dataSet.BarFillFac>1
        dataSet.BarFillFac=1;
    elseif dataSet.BarFillFac<0
        dataSet.BarFillFac=0;
    end
end
set(hObject,'String',mat2str(dataSet.BarFillFac));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function BarFillFacEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BarFillFacEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% Optimization Objectives

% --- Executes on button press in TorqueOptCheck.
function TorqueOptCheck_Callback(hObject, eventdata, handles)
% hObject    handle to TorqueOptCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TorqueOptCheck
dataSet = handles.dataSet;
dataSet.TorqueOptCheck = get(hObject,'Value');
% [bounds, objs, geo, per, mat] = data0(dataSet);
% dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
% handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

function MinExpTorEdit_Callback(hObject, eventdata, handles)
% Minimum expected torque (per.min_exp_torque)
% hObject    handle to MinExpTorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of MinExpTorEdit as text
%        str2double(get(hObject,'String')) returns contents of MinExpTorEdit as a double
dataSet = handles.dataSet;
dataSet.MinExpTorque = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function MinExpTorEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinExpTorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MaxExpeRippleTorEdit_Callback(hObject, eventdata, handles)
% Maximum expected torque ripple during optimization [p.u.](per.max_exp_ripple)
% hObject    handle to MaxExpeRippleTorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of MaxExpeRippleTorEdit as text
%        str2double(get(hObject,'String')) returns contents of MaxExpeRippleTorEdit as a double
dataSet = handles.dataSet;
dataSet.MaxRippleTorque = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function MaxExpeRippleTorEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxExpeRippleTorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in TorRipOptCheck.
function TorRipOptCheck_Callback(hObject, eventdata, handles)
% hObject    handle to TorRipOptCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TorRipOptCheck
dataSet = handles.dataSet;
dataSet.TorRipOptCheck = get(hObject,'Value');
% [bounds, objs, geo, per, mat] = data0(dataSet);
% dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
% handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)



function SlotSimulEdit_Callback(hObject, eventdata, handles)
% hObject    handle to SlotSimulEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SlotSimulEdit as text
%        str2double(get(hObject,'String')) returns contents of SlotSimulEdit as a double
dataSet = handles.dataSet;
Qs = str2double(get(hObject,'String'));
p = dataSet.NumOfPolePairs;
Q = round(dataSet.NumOfSlots*6*p);
yq = dataSet.PitchShortFac*dataSet.NumOfSlots*3;
path = pwd;
cd(fullfile (path,'koil'));
system(['koil_syre.exe',' ',num2str(Q),' ',num2str(p),' ',num2str(yq)]);
cd(path);
Windings = MatrixWin();
% check if Qs and ps is feasible (if ps is integer)
t2 = gcd(Q,2*p);
QsCalc = Q/t2;
psCalc = 2*p/t2;

ps = psCalc*Qs/QsCalc;

if not((ps-floor(ps))==0)
    disp('Error in the setting of Qs. Qs is reset to calculated standard value')
    ps=psCalc;
    Qs=QsCalc;
    set(handles.SlotSimulEdit,'String',int2str(Qs));
end
ps = round(ps);
dataSet.Qs = Qs;
%dataSet.ps = ps;
dataSet.WinMatr = Windings(:,1:floor(Qs)); % winding matrix, only Qs columns
columnName = cell(1,floor(Qs));
for i = 1 : floor(Qs)
    columnName{i} = ['Slot n° ',num2str(i)];
end
rowName{1} = 'Layer 1';
rowName{2} = 'Layer 2';
set(handles.WinTable,'rowname',rowName);
set(handles.WinTable,'columnname',columnName);
set(handles.WinTable,'data',dataSet.WinMatr(:,1:floor(Qs)));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function SlotSimulEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SlotSimulEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LossEvaluationCheck.
function LossEvaluationCheck_Callback(hObject, eventdata, handles)
% hObject    handle to LossEvaluationCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LossEvaluationCheck
dataSet = handles.dataSet;
dataSet.LossEvaluationCheck = get(hObject,'Value');
if dataSet.LossEvaluationCheck ==0
    set(handles.HysteresisLossFactorEdit,'Enable','off');
    set(handles.HysteresisFrequencyFactorEdit,'Enable','off');
    set(handles.HysteresisFluxDenEdit,'Enable','off');
    set(handles.EddyCurLossFactorEdit,'Enable','off');
    set(handles.MassDensityEdit,'Enable','off');
    set(handles.EvaluatedSpeedEdit,'Enable','off');
else
    [~, ~, ~, ~, mat] = data0(dataSet);
    set(handles.HysteresisLossFactorEdit,'Enable','on');
    set(handles.HysteresisFrequencyFactorEdit,'Enable','on');
    set(handles.HysteresisFluxDenEdit,'Enable','on');
    set(handles.EddyCurLossFactorEdit,'Enable','on');
    set(handles.MassDensityEdit,'Enable','on');
    set(handles.EvaluatedSpeedEdit,'Enable','on');
    dataSet.HysteresisLossFactor = mat.Rotor.kh;
    dataSet.HysteresisFrequencyFactor = mat.Rotor.alpha;
    dataSet.HysteresisFluxDenFactor = mat.Rotor.beta;
    dataSet.EddyCurLossFactorEdit = mat.Rotor.ke;
    dataSet.IronMassDen = mat.Rotor.kgm3;
    set(handles.HysteresisLossFactorEdit,'String',mat2str(dataSet.HysteresisLossFactor));
    set(handles.HysteresisFrequencyFactorEdit,'String',mat2str(dataSet.HysteresisFrequencyFactor));
    set(handles.HysteresisFluxDenEdit,'String',mat2str(dataSet.HysteresisFluxDenFactor));
    set(handles.EddyCurLossFactorEdit,'String',mat2str(dataSet.EddyCurLossFactorEdit));
    set(handles.MassDensityEdit,'String',mat2str(dataSet.IronMassDen));
%     Q = round(dataSet.NumOfSlots*6*dataSet.NumOfPolePairs);
%     t2 = gcd(Q,2*dataSet.NumOfPolePairs);
%     QsCalc = Q/t2;
%     dataSet.AngularSpanPP = max(60,20*QsCalc/dataSet.NumOfSlots);
%     set(handles.SpanEltPPEdit,'String',mat2str(dataSet.AngularSpanPP));
%     dataSet.NumOfRotPosPP = ceil(dataSet.AngularSpanPP/60)*6;
%     set(handles.NumOfRotorPosiPPEdit,'String',mat2str(dataSet.NumOfRotPosPP));
    handles.dataSet = dataSet;
    dataSet = DrawPushMachine(handles,dataSet.currentfilename,dataSet.currentpathname);
end
handles.dataSet = dataSet;
guidata(hObject,handles)



function HysteresisLossFactorEdit_Callback(hObject, eventdata, handles)
% hObject    handle to HysteresisLossFactorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HysteresisLossFactorEdit as text
%        str2double(get(hObject,'String')) returns contents of HysteresisLossFactorEdit as a double
dataSet = handles.dataSet;
dataSet.HysteresisLossFactor = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function HysteresisLossFactorEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HysteresisLossFactorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function HysteresisFrequencyFactorEdit_Callback(hObject, eventdata, handles)
% hObject    handle to HysteresisFrequencyFactorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HysteresisFrequencyFactorEdit as text
%        str2double(get(hObject,'String')) returns contents of HysteresisFrequencyFactorEdit as a double
dataSet = handles.dataSet;
dataSet.HysteresisFrequencyFactor = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function HysteresisFrequencyFactorEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HysteresisFrequencyFactorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function HysteresisFluxDenEdit_Callback(hObject, eventdata, handles)
% hObject    handle to HysteresisFluxDenEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HysteresisFluxDenEdit as text
%        str2double(get(hObject,'String')) returns contents of HysteresisFluxDenEdit as a double
dataSet = handles.dataSet;
dataSet.HysteresisFluxDenFactor = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function HysteresisFluxDenEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HysteresisFluxDenEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EddyCurLossFactorEdit_Callback(hObject, eventdata, handles)
% hObject    handle to EddyCurLossFactorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EddyCurLossFactorEdit as text
%        str2double(get(hObject,'String')) returns contents of EddyCurLossFactorEdit as a double
dataSet = handles.dataSet;
dataSet.EddyCurLossFactor = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function EddyCurLossFactorEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EddyCurLossFactorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MassDensityEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MassDensityEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MassDensityEdit as text
%        str2double(get(hObject,'String')) returns contents of MassDensityEdit as a double
dataSet = handles.dataSet;
dataSet.IronMassDen = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function MassDensityEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MassDensityEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EvaluatedSpeedEdit_Callback(hObject, eventdata, handles)
% hObject    handle to EvaluatedSpeedEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EvaluatedSpeedEdit as text
%        str2double(get(hObject,'String')) returns contents of EvaluatedSpeedEdit as a double
dataSet = handles.dataSet;
dataSet.EvalSpeed = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function EvaluatedSpeedEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EvaluatedSpeedEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MassCuEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MassCuEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MassCuEdit as text
%        str2double(get(hObject,'String')) returns contents of MassCuEdit as a double

dataSet = handles.dataSet;
dataSet.MaxCuMass = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function MassCuEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MassCuEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MassCuOptCheck.
function MassCuOptCheck_Callback(hObject, eventdata, handles)
% hObject    handle to MassCuOptCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MassCuOptCheck

dataSet = handles.dataSet;
dataSet.MassCuOptCheck = get(hObject,'Value');
% [bounds, objs, geo, per, mat] = data0(dataSet);
% dataSet.RQ = buildDefaultRQ(bounds);
handles.dataSet = dataSet;
% handles = DrawPush_Callback(hObject, eventdata, handles);
guidata(hObject,handles)


% --- Executes on button press in SlotLayerPosCheck.
function SlotLayerPosCheck_Callback(hObject, eventdata, handles)
% hObject    handle to SlotLayerPosCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SlotLayerPosCheck

dataSet = handles.dataSet;
dataSet.SlotLayerPosCheck = get(hObject,'Value');
% [bounds, objs, geo, per, mat] = data0(dataSet);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject,eventdata,handles);
guidata(hObject,handles)


% --- Executes on button press in RadRibCheck.
function RadRibCheck_Callback(hObject, eventdata, handles)
% hObject    handle to RadRibCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RadRibCheck

dataSet = handles.dataSet;
dataSet.RadRibCheck = get(hObject,'Value');
% [bounds, objs, geo, per, mat] = data0(dataSet);
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject,eventdata,handles);
dataSet.RadRibEdit = eval(get(handles.RadRibEdit,'String'));
handles.dataSet = dataSet;
SetParameters(handles,dataSet)
guidata(hObject,handles)



function RadRibEdit_Callback(hObject, eventdata, handles)
% hObject    handle to RadRibEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RadRibEdit as text
%        str2double(get(hObject,'String')) returns contents of RadRibEdit as a double

dataSet = handles.dataSet;
dataSet.RadRibEdit = eval(get(hObject,'String'));
handles.dataSet = dataSet;
handles = DrawPush_Callback(hObject,eventdata,handles);
SetParameters(handles,dataSet)
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function RadRibEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RadRibEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in IronLibraryPush.
function IronLibraryPush_Callback(hObject, eventdata, handles)
% hObject    handle to IronLibraryPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mat = material_properties_iron(0);
material = char(mat.MatList{1});
for ii=2:length(mat.MatList)
    material = [material ', ' char(mat.MatList{ii})];
end

set(handles.MaterialText,'Style','Edit');
set(handles.MaterialText,'String',material);


% --- Executes on button press in ConductorLibraryPush.
function ConductorLibraryPush_Callback(hObject, eventdata, handles)
% hObject    handle to ConductorLibraryPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mat = material_properties_conductor(0);
material = char(mat.MatList{1});
for ii=2:length(mat.MatList)
    material = [material ', ' char(mat.MatList{ii})];
end

set(handles.MaterialText,'Style','Edit');
set(handles.MaterialText,'String',material);

% --- Executes on button press in BarrierLibraryPush.
function BarrierLibraryPush_Callback(hObject, eventdata, handles)
% hObject    handle to BarrierLibraryPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mat = material_properties_layer(0);
material = char(mat.MatList{1});
for ii=2:length(mat.MatList)
    material = [material ', ' char(mat.MatList{ii})];
end

set(handles.MaterialText,'Style','Edit');
set(handles.MaterialText,'String',material);


% --- Executes on button press in ViewPropPush.
function ViewPropPush_Callback(hObject, eventdata, handles)
% hObject    handle to ViewPropPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

answer = inputdlg('Material Name','Properties viewer',1);
MatName = char(answer);
mat = material_properties_iron(MatName);
if isfield(mat,'kgm3')
    plot_iron_prop(mat);
else
    clear mat
    mat = material_properties_conductor(MatName);
    if isfield(mat,'kgm3')
        plot_conductor_prop(mat);
    else
        clear mat
        mat = material_properties_layer(MatName);
        if isfield(mat,'kgm3')
            plot_layer_prop(mat);
        else
            msgbox('Select a correct material name to see the properties','Properties viewer');
        end
    end
end



function BfeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to BfeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BfeEdit as text
%        str2double(get(hObject,'String')) returns contents of BfeEdit as a double
dataSet = handles.dataSet;
dataSet.Bfe = str2num(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function BfeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BfeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ktEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ktEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ktEdit as text
%        str2double(get(hObject,'String')) returns contents of ktEdit as a double
dataSet = handles.dataSet;
dataSet.kt = str2num(get(hObject,'String'));
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ktEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ktEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
