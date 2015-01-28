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
%      GUI_SYRE, by itself, creates a new GUI_SYRE or raises the existing
%      singleton*.
%
%      H = GUI_SYRE returns the handle to a new GUI_SYRE or the handle to
%      the existing singleton*.
%
%      GUI_SYRE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SYRE.M with the given input arguments.
%
%      GUI_SYRE('Property','Value',...) creates a new GUI_SYRE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Syre_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Syre_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Syre

% Last Modified by GUIDE v2.5 14-Oct-2014 12:29:31

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
% addpath('./geo');

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
handles.left_right_corner = handles.pan1pos(2) + handles.pan1pos(4);
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

%% === FLAG ===============================================================
handles.MatrixWinFlag = 1;
handles.UpdateData0 = 0;
handles.Opti = 0;
%% ========================================================================


% Update handles structure
guidata(hObject, handles);
axes(handles.axes4);
SyreImg = imread('syre.png'); 
image(SyreImg);
axis off;
axes(handles.axes5);
set(handles.axes5,'YTickLabel',[]);
set(handles.axes5,'XTickLabel',[]);
box on;

% UIWAIT makes GUI_Syre wait for user response (see UIRESUME)
% uiwait(handles.figure1);
[bounds, geo, per] = data0();

dataSet.NumOfPolePairs = geo.p; % number of pole pairs
dataSet.AirGapThickness = geo.g; % airgap thickness
dataSet.StatorOuterRadius = geo.R; % stator outer radius
dataSet.AirGapRadius = geo.r; % machine airgap radius
dataSet.ShaftRadius = geo.Ar; % shaft radius
dataSet.StackLength = geo.l; % stack length
dataSet.TypeOfRotor = geo.RotType; % type of rotor
dataSet.NumOfSlots = geo.q; % number of slots
dataSet.ToothLength = geo.lt; % tooth length
dataSet.StatorSlotOpen = geo.acs; % stator slot open in [p.u.]
dataSet.ToothWidth = geo.wt; % Bgap/Btooth (determines tooth width^-1, yoke width^-1)
dataSet.ToothTangDepth = geo.ttd; % tooth tang depth [mm]
dataSet.ToothTangAngle = geo.tta; %  tooth tang angle (mech degree)
dataSet.FilletCorner = geo.SFR; % fillet at the back corner of the slot [mm]
dataSet.SlotMaterial = geo.BLKLABELSmaterials{3}; % slot material
dataSet.StatorMaterial = geo.BLKLABELSmaterials{4}; % stator material
dataSet.RotorMaterial = geo.BLKLABELSmaterials{5}; % rotor material
dataSet.FluxBarrierMaterial = geo.BLKLABELSmaterials{6}; % flux barrier material
dataSet.ShaftMaterial = geo.BLKLABELSmaterials{7}; % shaft material
dataSet.RotorCondMaterial = geo.BLKLABELSmaterials{8}; % rotor conductor material
dataSet.SlotFillFactor = geo.kcu; % slot fill factor
dataSet.PitchShortFac = geo.kracc; % pitch short factor
dataSet.TurnsInSeries= geo.Ns; % turns in series
dataSet.AdmiJouleLosses = per.Loss; % admitted Joule Losses [W]
dataSet.CopperTemp = per.tempcu; % copper temperature [C]
dataSet.CurrOverLoad = per.overload; % Current overload
dataSet.NumOfLayers = geo.nlay; % Number of Layers (bars)
dataSet.OverSpeed = geo.nmax; % overspeed in [rpm]
dataSet.Br = geo.Br; % Br
dataSet.Hc = geo.Hc; % Hc
dataSet.MinExpTorque = per.min_exp_torque; % Minimum expected torque
dataSet.MaxRippleTorque = per.max_exp_ripple; % Maximum ripple torque
dataSet.MinMechTol = geo.pont0; % minimum mechanical tolerance [mm]
dataSet.SimPoMOOA = geo.nsim_MOOA; % simulated positions (6-1)
dataSet.RotPoMOOA = geo.delta_sim_MOOA; % rotor position span [elt degrees]
dataSet.SimPoFine = geo.nsim_singt; % simulated positions (16-1)
dataSet.RotPoFine = geo.delta_sim_singt; % rotor position span [elt degrees]
dataSet.Mesh = geo.K_mesh; % Mesh
dataSet.Mesh_MOOA = geo.K_mesh_MOOA; % Mesh MOOA
dataSet.RQnames = geo.RQnames;

%% ==== BOUNDS ============================================================
dataSet.Alpha1Bou = [22.5 45]/90;
dataSet.DeltaAlphaBou = round([0.5/geo.nlay 0.5] * 100)/100; % other angles [p.u.]
dataSet.hcBou = [0.2 1];               % barrier ticknesses [p.u.]
dataSet.DfeBou = [-0.75 0.75];         % barrier offset [p.u.]
dataSet.BrBou = [0.2 0.8]; 
dataSet.GapBou = rrtd(dataSet.AirGapThickness*[0.7  1.5],-1);
dataSet.GapRadiusBou = rrtd(dataSet.AirGapRadius*[0.8  1.2],-1);
dataSet.ToothWiBou = rrtd(dataSet.ToothWidth*[0.75  1.3],-1);
dataSet.ToothLeBou = rrtd(dataSet.ToothLength*[0.8 1.2],-1);
dataSet.PhaseAngleCurrBou = bounds(end,:);          % phase angle of the current vector
dataSet.Dalpha1BouCheck = 1;
dataSet.DalphaBouCheck = 1;
dataSet.hcBouCheck = 1;
if (strcmp(geo.RotType,'Fluid') || strcmp(geo.RotType,'Seg'))
    dataSet.DxBouCheck = 1;
else
    dataSet.DxBouCheck = 0;
end
dataSet.GammaBouCheck = 1;
dataSet.GapBouCheck = 0;
dataSet.BrBouCheck = 0;
dataSet.AirgapRadiusBouCheck = 0;
dataSet.ToothWidthBouCheck = 0;
dataSet.ToothLengthBouCheck = 0;

%% ==== OPTIMIZATION ======================================================
dataSet.MaxGen = 3;
dataSet.XPop = 4;

%% ========================================================================
%% ==== For Post Processing ===============================================

dataSet.CurrLoPP = 1; % [p.u.] current load post processing
dataSet.GammaPP = 60; % [deg] current load post processing
dataSet.BrPP = 0; % [T] Br
dataSet.NumGrid = 5; % number of points in [0 Imax] for the single machine post-processing
% dataSet.NumInter = per.n_interp;  % number of points in [0 Imax] for data interpolation
dataSet.NumOfRotPosPP = 20;
dataSet.AngularSpanPP = 60;
%% ==== RQ per plot =======================================================
% data = buildDefaultRQ(dataSet,bounds);
% set(handles.RQPlotEdit,'String',mat2str(data));
% dataSet.RQ = data;

%% ========================================================================
dataSet.XFEMMOpt = 'N';
dataSet.XFEMMPPMot = 'N';
% if exist('bounds_Delta_X')
%     dataSet.DeltaXBou = bounds_Delta_X; % p.u. displacement of the 1st layer (Seg only)
% else
%     dataSet.DeltaXBou = [1 1];
% end
dataSet.RMVTmp = geo.RemoveTMPfile; %  for removing the motor folders in tmp
%% ==== Matrix of winding =================================================
p = dataSet.NumOfPolePairs;
Q = round(dataSet.NumOfSlots*6*p);
yq = dataSet.PitchShortFac*dataSet.NumOfSlots*3;
path = pwd;
cd(fullfile (path,'koil'));
system(['koil_syre.exe',' ',num2str(Q),' ',num2str(p),' ',num2str(yq)]);
cd(path);
Windings = MatrixWin();
%dataSet.WinMatr = Windings; % winding matrix
t = gcd(round(dataSet.NumOfSlots*6*dataSet.NumOfPolePairs),dataSet.NumOfPolePairs);  % periodicity
if ((6*t/Q)>1)
    Qs = Q/t;   % periodic machine
else
    Qs = Q/2/t; % anti-periodic machine
end
dataSet.WinMatr = Windings(:,1:floor(Qs)); % winding matrix, only Qs columns
%% ==== RQ per plot =======================================================
data = buildDefaultRQ(dataSet,bounds);
% set(handles.RQPlotEdit,'String',mat2str(data));
dataSet.RQ = data;
dataSet.HCpu = data((dataSet.NumOfLayers+1):end-1);
dataSet.ALPHApu = data(1:dataSet.NumOfLayers);

%% ========================================================================
handles.dataSet = dataSet; % data structure
SetParameters(handles,dataSet) % aux. function for set the values in the edit boxes
%% ======= Matrix of winding visual =======================================
p = dataSet.NumOfPolePairs;
Q = dataSet.NumOfSlots*6*p;
%dataSet.WinMatr = Windings; % winding matrix
t = gcd(round(dataSet.NumOfSlots*6*dataSet.NumOfPolePairs),dataSet.NumOfPolePairs);  % periodicity
if ((6*t/Q)>1)
    Qs = Q/t;   % periodic machine
else
    Qs = Q/2/t; % anti-periodic machine
end
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
path = pwd
cd(fullfile (path,'koil'));
system(['koil_syre.exe',' ',num2str(Q),' ',num2str(p),' ',num2str(yq)]);
cd(path);
Windings = MatrixWin();
%dataSet.WinMatr = Windings; % winding matrix
t = gcd(round(dataSet.NumOfSlots*6*dataSet.NumOfPolePairs),dataSet.NumOfPolePairs);  % periodicity
if ((6*t/Q)>1)
    Qs = Q/t;   % periodic machine
else
    Qs = Q/2/t; % anti-periodic machine
end
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
DrawPush_Callback(hObject, eventdata, handles);
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
set(handles.CopperTempEdit,'String',num2str(dataIn.CopperTemp));
% set(handles.DCLinkVolEdit,'String',num2str(dataIn.DCVoltage));
set(handles.CurrentOverLoadEdit,'String',num2str(dataIn.CurrOverLoad));
set(handles.NumberOfLayersEdit,'String',num2str(dataIn.NumOfLayers));
set(handles.OverSpeedEdit,'String',num2str(dataIn.OverSpeed));
set(handles.BrPMEdit,'String',num2str(dataIn.Br));
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
set(handles.BrPPEdit,'String',mat2str(dataIn.BrPP));
set(handles.CurrLoPPEdit,'String',mat2str(dataIn.CurrLoPP));
set(handles.GammaPPEdit,'String',mat2str(dataIn.GammaPP));
% set(handles.RQPlotEdit,'String',mat2str(dataIn.RQ));
set(handles.AlphapuEdit,'String',mat2str(dataIn.ALPHApu));
set(handles.hcpuEdit,'String',mat2str(dataIn.HCpu));
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
%% === Values for the edit text Plot ======================================
flag_plot = 'Y';
h = handles.axes5;
[hc,dalpha] = Plot_Machine(h,dataIn,flag_plot);
view = round(100*[dalpha hc])/100;
% set(handles.ViewHcDaEdit,'String',mat2str(view));
set(handles.AlphapuEdit,'String',mat2str(dataIn.RQ(1:dataIn.NumOfLayers)));
set(handles.hcpuEdit,'String',mat2str(dataIn.RQ((dataIn.NumOfLayers+1):end-1)));
set(handles.AlphadegreeEdit,'String',mat2str(view(1:dataIn.NumOfLayers)));
set(handles.hcmmEdit,'String',mat2str(view((dataIn.NumOfLayers+1):end)));



function GapThiEdit_Callback(hObject, eventdata, handles)
% hObject    handle to GapThiEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of GapThiEdit as text
%        str2double(get(hObject,'String')) returns contents of GapThiEdit as a double
dataSet = handles.dataSet;
dataSet.AirGapThickness = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
DrawPush_Callback(hObject, eventdata, handles);
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
% hObject    handle to StatorOuterRadEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of StatorOuterRadEdit as text
%        str2double(get(hObject,'String')) returns contents of StatorOuterRadEdit as a double
dataSet = handles.dataSet;
dataSet.StatorOuterRadius = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
DrawPush_Callback(hObject, eventdata, handles);
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
% hObject    handle to AirGapRadiusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of AirGapRadiusEdit as text
%        str2double(get(hObject,'String')) returns contents of AirGapRadiusEdit as a double
dataSet = handles.dataSet;
dataSet.AirGapRadius = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
DrawPush_Callback(hObject, eventdata, handles);
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
% hObject    handle to ShaftRadEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of ShaftRadEdit as text
%        str2double(get(hObject,'String')) returns contents of ShaftRadEdit as a double
dataSet = handles.dataSet;
dataSet.ShaftRadius = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
DrawPush_Callback(hObject, eventdata, handles);
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
% hObject    handle to StackLenghtEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of StackLenghtEdit as text
%        str2double(get(hObject,'String')) returns contents of StackLenghtEdit as a double
dataSet = handles.dataSet;
dataSet.StackLength = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
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
% hObject    handle to TypeOfRotorList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns TypeOfRotorList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TypeOfRotorList
contents = get(hObject,'String');
s = contents{get(hObject,'Value')};
dataSet = handles.dataSet;
dataSet.TypeOfRotor = s;
[bounds, geo, per] = data0(dataSet);
data = buildDefaultRQ(dataSet,bounds);
set(handles.hcpuEdit,'String',mat2str(data((dataSet.NumOfLayers+1):end-1)));
set(handles.AlphapuEdit,'String',mat2str(data(1:dataSet.NumOfLayers)));
[truefalse, index] = ismember('dx', geo.RQnames);
if truefalse 
    dataSet.DfeBou = bounds(index,:);            % barrier offset [p.u.]
end
dataSet.RQ = data;
flag_plot = 'Y';
h = handles.axes5;
[hc,dalpha] = Plot_Machine(h,dataSet,flag_plot);
view = round(100*[dalpha hc])/100;
set(handles.AlphadegreeEdit,'String',mat2str(view(1:dataSet.NumOfLayers)));
set(handles.hcmmEdit,'String',mat2str(view((dataSet.NumOfLayers+1):end)));
set(handles.DfeBouEdit,'String',mat2str(dataSet.DfeBou));
handles.dataSet = dataSet;
DrawPush_Callback(hObject, eventdata, handles);
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
t = gcd(round(dataSet.NumOfSlots*6*dataSet.NumOfPolePairs),dataSet.NumOfPolePairs);  % periodicity
if ((6*t/Q)>1)
    Qs = Q/t;   % periodic machine
else
    Qs = Q/2/t; % anti-periodic machine
end
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
DrawPush_Callback(hObject, eventdata, handles);
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
% hObject    handle to ToothLengEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of ToothLengEdit as text
%        str2double(get(hObject,'String')) returns contents of ToothLengEdit as a double
dataSet = handles.dataSet;
dataSet.ToothLength = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
DrawPush_Callback(hObject, eventdata, handles);
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
% hObject    handle to StatorSlotOpeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of StatorSlotOpeEdit as text
%        str2double(get(hObject,'String')) returns contents of StatorSlotOpeEdit as a double
dataSet = handles.dataSet;
dataSet.StatorSlotOpen = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
DrawPush_Callback(hObject, eventdata, handles);
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
% hObject    handle to ToothWidthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of ToothWidthEdit as text
%        str2double(get(hObject,'String')) returns contents of ToothWidthEdit as a double
dataSet = handles.dataSet;
dataSet.ToothWidth = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
DrawPush_Callback(hObject, eventdata, handles);
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
% hObject    handle to ToothTanDepEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of ToothTanDepEdit as text
%        str2double(get(hObject,'String')) returns contents of ToothTanDepEdit as a double
dataSet = handles.dataSet;
dataSet.ToothTangDepth = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
DrawPush_Callback(hObject, eventdata, handles);
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
% hObject    handle to ToothTangAngleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of ToothTangAngleEdit as text
%        str2double(get(hObject,'String')) returns contents of ToothTangAngleEdit as a double
dataSet = handles.dataSet;
dataSet.ToothTangAngle = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
DrawPush_Callback(hObject, eventdata, handles);
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
% hObject    handle to FillCorSlotEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of FillCorSlotEdit as text
%        str2double(get(hObject,'String')) returns contents of FillCorSlotEdit as a double
dataSet = handles.dataSet;
dataSet.FilletCorner = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
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
% hObject    handle to SlotFillFacEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of SlotFillFacEdit as text
%        str2double(get(hObject,'String')) returns contents of SlotFillFacEdit as a double
dataSet = handles.dataSet;
dataSet.SlotFillFactor = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
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
t = gcd(round(dataSet.NumOfSlots*6*dataSet.NumOfPolePairs),dataSet.NumOfPolePairs);  % periodicity
if ((6*t/Q)>1)
    Qs = Q/t;   % periodic machine
else
    Qs = Q/2/t; % anti-periodic machine
end
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
% hObject    handle to TurnsSeriesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of TurnsSeriesEdit as text
%        str2double(get(hObject,'String')) returns contents of TurnsSeriesEdit as a double
dataSet = handles.dataSet;
dataSet.TurnsInSeries = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
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
% hObject    handle to JouleLossesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of JouleLossesEdit as text
%        str2double(get(hObject,'String')) returns contents of JouleLossesEdit as a double
dataSet = handles.dataSet;
dataSet.AdmiJouleLosses = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
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
% hObject    handle to CopperTempEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of CopperTempEdit as text
%        str2double(get(hObject,'String')) returns contents of CopperTempEdit as a double
dataSet = handles.dataSet;
dataSet.CopperTemp = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
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
% hObject    handle to NumberOfLayersEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of NumberOfLayersEdit as text
%        str2double(get(hObject,'String')) returns contents of NumberOfLayersEdit as a double
dataSet = handles.dataSet;
dataSet.NumOfLayers = str2double(get(hObject,'String'));
[bounds, geo, per] = data0(dataSet);
dataSet.RQnames = geo.RQnames;
data = buildDefaultRQ(dataSet,bounds);
set(handles.hcpuEdit,'String',mat2str(data((dataSet.NumOfLayers+1):end-1)));
set(handles.AlphapuEdit,'String',mat2str(data(1:dataSet.NumOfLayers)));
dataSet.RQ = data;
flag_plot = 'Y';
h = handles.axes5;
[hc,dalpha] = Plot_Machine(h,dataSet,flag_plot);
view = round(100*[dalpha hc])/100;
set(handles.hcmmEdit,'String',mat2str(view((dataSet.NumOfLayers+1):end)));
set(handles.AlphadegreeEdit,'String',mat2str(view(1:dataSet.NumOfLayers)));
set(handles.DfeBouEdit,'String',mat2str(dataSet.DfeBou));
handles.dataSet = dataSet;
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
% hObject    handle to OverSpeedEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of OverSpeedEdit as text
%        str2double(get(hObject,'String')) returns contents of OverSpeedEdit as a double
dataSet = handles.dataSet;
dataSet.OverSpeed = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
DrawPush_Callback(hObject, eventdata, handles);
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
% hObject    handle to BrPMEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of BrPMEdit as text
%        str2double(get(hObject,'String')) returns contents of BrPMEdit as a double
dataSet = handles.dataSet;
dataSet.Br = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
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

function MinExpTorEdit_Callback(hObject, eventdata, handles)
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

function MecTolerEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MecTolerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of MecTolerEdit as text
%        str2double(get(hObject,'String')) returns contents of MecTolerEdit as a double
dataSet = handles.dataSet;
dataSet.MinMechTol = str2double(get(hObject,'String'));
handles.dataSet = dataSet;
DrawPush_Callback(hObject, eventdata, handles);
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
% hObject    handle to FluxBarMatEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of FluxBarMatEdit as text
%        str2double(get(hObject,'String')) returns contents of FluxBarMatEdit as a double
dataSet = handles.dataSet;
dataSet.FluxBarrierMaterial = get(hObject,'String');
handles.dataSet = dataSet;
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
function LibraryMaterPush_Callback(hObject, eventdata, handles)
% hObject    handle to LibraryMaterPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileID = fopen('empty_case.fem');
C = textscan(fileID,'%s');
fclose(fileID);
n = length(C{1});
A = zeros(100,2);
j = 0;
for i = 3 : 1 : n
    if strcmp(C{1}(i-2),'<BlockName>') && strcmp(C{1}(i-1),'=')
        j = j + 1;
        A(j,1) = i; %matrice degli indici dei materiali
        while ~strcmp(C{1}(i+1),'<Mu_x>')
            i = i + 1;
        end
        A(j,2) = i;
    end
end
[m,~] = size(A);
for i = 1 : m
    if A(i,:) == 0
        m = i - 1;
        break
    end
end
material = {};
for i = 1 : m
    if A(i,1) == A(i,2)
        material = strjoin([material,C{1}(A(i))]);
    else
        diff = A(i,2) - A(i,1);
        index = A(i,1);
        name = {};
        for j = 0 : diff
            name = strcat([name,C{1}(index + j)]);
        end
        material = strjoin([material,name]);
    end
end

set(handles.MaterialText,'Style','Edit');
set(handles.MaterialText,'Max',m);
set(handles.MaterialText,'String',material);
guidata(hObject,handles)

% --- Executes on button press in OptimizePush.
function OptimizePush_Callback(hObject, eventdata, handles)
% hObject    handle to OptimizePush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dataSet = handles.dataSet;
figure()
[bounds, geo, per] = data0(dataSet);
dat.geo0 = geo;
dat.per = per;
%%%%%%%%%% FEMM fitness handle %%%%%%%%%%%%%%%%%%%%%%%%%%
eval_type = 'MO_OA';
if strcmp(dataSet.XFEMMOpt,'Y')
    FitnessFunction = @(x)FEMMfitnessX(x,geo,per,eval_type);
else
    FitnessFunction = @(x)FEMMfitness(x,geo,per,eval_type);
end
dat.CostProblem = FitnessFunction;           % Cost function instance
%% ========================================================================
NOBJ = 2;
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
OUT = MODE2(dat);
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
function SavePushTool_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to SavePushTool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dataSet = handles.dataSet;
uisave('dataSet','Project.mat');
guidata(hObject,handles);

% --------------------------------------------------------------------
function LoadPushTool_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to LoadPushTool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.mat');
load([PathName FileName]);
SetParameters(handles,dataSet);
handles.dataSet = dataSet;
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
set(handles.RotorGeometryPanel,'Visible','off');

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
set(handles.RotorGeometryPanel,'Visible','off');

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
                [bounds, geo] = data0();
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
% hObject    handle to WinTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in DrawPush.
function DrawPush_Callback(hObject, eventdata, handles)
% hObject    handle to DrawPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
flag_plot = 'Y';
h = handles.axes5;
dataSet = handles.dataSet;
[hc,dalpha] = Plot_Machine(h,dataSet,flag_plot);
view = round(100*[dalpha hc])/100;
set(handles.ViewHcDaEdit,'String',mat2str(view));
guidata(hObject,handles)

function MeshEdit_Callback(hObject, eventdata, handles)
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
handles.dataSet = dataSet;
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
handles.dataSet = dataSet;
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
handles.dataSet = dataSet;
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
handles.dataSet = dataSet;
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
% hObject    handle to BrPPEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of BrPPEdit as text
%        str2double(get(hObject,'String')) returns contents of BrPPEdit as a double

dataSet = handles.dataSet;
try
    dataSet.BrPP = str2num((get(hObject,'String')));
catch
    dataSet.BrPP = str2double((get(hObject,'String')));
end
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
% hObject    handle to StartPProPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataSet = handles.dataSet;
if strcmp(dataSet.XFEMMPPMot,'Y')
    post_proc_single_motorX(dataSet.CurrLoPP,dataSet.GammaPP,dataSet.BrPP,dataSet.NumOfRotPosPP,dataSet.AngularSpanPP,dataSet.NumGrid);
else
    post_proc_single_motor(dataSet.CurrLoPP,dataSet.GammaPP,dataSet.BrPP,dataSet.NumOfRotPosPP,dataSet.AngularSpanPP,dataSet.NumGrid);
end
guidata(hObject,handles)


% --- Executes on button press in PostProcXFEMMCheck.
function PostProcXFEMMCheck_Callback(hObject, eventdata, handles)
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
% hObject    handle to MatrixOfData0Push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of MatrixOfData0Push
dataSet = handles.dataSet;
if exist([cd,'\tmp\flag.mat']) > 0
    flag = 0;
    save([cd,'\tmp\flag.mat'],'flag');
    [bounds, geo] = data0();    
    flag = 1;
    save([cd,'\tmp\flag.mat'],'flag');
else
    [bounds, geo] = data0();  
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
% hObject    handle to SaveMachinePush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DrawPushMachine;
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
% hObject    handle to AlphapuEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of AlphapuEdit as text
%        str2double(get(hObject,'String')) returns contents of AlphapuEdit as a double
dataSet = handles.dataSet;
dataSet.ALPHApu = str2num(get(hObject,'String'));
dataSet.HCpu = str2num(get(handles.hcpuEdit,'String'));
dataSet.RQ = [dataSet.ALPHApu dataSet.HCpu];
handles.dataSet = dataSet;
flag_plot = 'Y';
h = handles.axes5;
[hc,dalpha] = Plot_Machine(h,dataSet,flag_plot);
view = round(100*[dalpha hc])/100;
set(handles.AlphadegreeEdit,'String',mat2str(view(1:dataSet.NumOfLayers)));
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
% hObject    handle to AlphadegreeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of AlphadegreeEdit as text
%        str2double(get(hObject,'String')) returns contents of AlphadegreeEdit as a double

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
% hObject    handle to hcpuEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of hcpuEdit as text
%        str2double(get(hObject,'String')) returns contents of hcpuEdit as a double
dataSet = handles.dataSet;
dataSet.HCpu = str2num(get(hObject,'String'));
dataSet.ALPHApu = str2num(get(handles.AlphapuEdit,'String'));
dataSet.RQ = [dataSet.ALPHApu dataSet.HCpu];
handles.dataSet = dataSet;
h = handles.axes5;
flag_plot = 'Y';
[hc,dalpha] = Plot_Machine(h,dataSet,flag_plot);
view = round(100*[dalpha hc])/100;
set(handles.hcmmEdit,'String',mat2str(view((dataSet.NumOfLayers+1):end)));
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
% hObject    handle to hcmmEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of hcmmEdit as text
%        str2double(get(hObject,'String')) returns contents of hcmmEdit as a double


% --- Executes during object creation, after setting all properties.
function hcmmEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hcmmEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BrBouEdit_Callback(hObject, eventdata, handles)
% hObject    handle to BrBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BrBouEdit as text
%        str2double(get(hObject,'String')) returns contents of BrBouEdit as a double
dataSet = handles.dataSet;
dataSet.BrBou = str2num(get(hObject,'String'));
handles.dataSet = dataSet;
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
% hObject    handle to BrBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of BrBouCheck
dataSet = handles.dataSet;
dataSet.BrBouCheck = get(hObject,'Value');
handles.dataSet = dataSet;
guidata(hObject,handles)




function AirgapRadiusBouEdit_Callback(hObject, eventdata, handles)
% hObject    handle to AirgapRadiusBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AirgapRadiusBouEdit as text
%        str2double(get(hObject,'String')) returns contents of AirgapRadiusBouEdit as a double
dataSet = handles.dataSet;
dataSet.GapRadiusBou = str2num(get(hObject,'String'));
handles.dataSet = dataSet;
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
% hObject    handle to AirgapRadiusBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AirgapRadiusBouCheck
dataSet = handles.dataSet;
dataSet.AirgapRadiusBouCheck = get(hObject,'Value');
handles.dataSet = dataSet;
guidata(hObject,handles)



function ToothWidthBouEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ToothWidthBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ToothWidthBouEdit as text
%        str2double(get(hObject,'String')) returns contents of ToothWidthBouEdit as a double
dataSet = handles.dataSet;
dataSet.ToothWiBou = str2num(get(hObject,'String'));
handles.dataSet = dataSet;
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
% hObject    handle to ToothWidthBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ToothWidthBouCheck
dataSet = handles.dataSet;
dataSet.ToothWidthBouCheck = get(hObject,'Value');
handles.dataSet = dataSet;
guidata(hObject,handles)


function ToothLenBouEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ToothLenBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ToothLenBouEdit as text
%        str2double(get(hObject,'String')) returns contents of ToothLenBouEdit as a double
dataSet = handles.dataSet;
dataSet.ToothLeBou = str2num(get(hObject,'String'));
handles.dataSet = dataSet;
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
% hObject    handle to ToothLengthBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ToothLengthBouCheck
dataSet = handles.dataSet;
dataSet.ToothLengthBouCheck = get(hObject,'Value');
handles.dataSet = dataSet;
guidata(hObject,handles)


function GapBouEdit_Callback(hObject, eventdata, handles)
% hObject    handle to GapBouEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GapBouEdit as text
%        str2double(get(hObject,'String')) returns contents of GapBouEdit as a double


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
% hObject    handle to Dalpha1BouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Dalpha1BouCheck
dataSet = handles.dataSet;
dataSet.Dalpha1BouCheck = get(hObject,'Value');
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes on button press in DalphaBouCheck.
function DalphaBouCheck_Callback(hObject, eventdata, handles)
% hObject    handle to DalphaBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DalphaBouCheck
dataSet = handles.dataSet;
dataSet.DalphaBouCheck = get(hObject,'Value');
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes on button press in hcBouCheck.
function hcBouCheck_Callback(hObject, eventdata, handles)
% hObject    handle to hcBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hcBouCheck
dataSet = handles.dataSet;
dataSet.hcBouCheck = get(hObject,'Value');
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes on button press in DxBouCheck.
function DxBouCheck_Callback(hObject, eventdata, handles)
% hObject    handle to DxBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DxBouCheck
dataSet = handles.dataSet;
dataSet.DxBouCheck = get(hObject,'Value');
handles.dataSet = dataSet;
guidata(hObject,handles)


% --- Executes on button press in GammaBouCheck.
function GammaBouCheck_Callback(hObject, eventdata, handles)
% hObject    handle to GammaBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GammaBouCheck
dataSet = handles.dataSet;
dataSet.GammaBouCheck = get(hObject,'Value');
handles.dataSet = dataSet;
guidata(hObject,handles)

% --- Executes on button press in GapBouCheck.
function GapBouCheck_Callback(hObject, eventdata, handles)
% hObject    handle to GapBouCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of GapBouCheck
dataSet = handles.dataSet;
dataSet.GapBouCheck = get(hObject,'Value');
handles.dataSet = dataSet;
guidata(hObject,handles)
