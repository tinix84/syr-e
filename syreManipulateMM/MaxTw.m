% Copyright 2015
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%
%%%%%%%%%
% MaxTw.m
% Control trajectories for loss minimization in the T-w plane
% Eval loss and effy map accordingly

% Input:
% 1) fdfq_idiq_n256.mat    flux map (id,iq) 
% 2) fdfq_idiq_loss.mat    iron loss map on id iq (optional)
% 3) ReadParameters.m      general parameters 
% [(1) and (2) can be a single file, load it twice when prompted for]

% Notes:
% - Voltage limit is respected (parameter d.Vbus in ReadParamentrs.m)
% - Current limit is not imposed but visible in the current amplitude map
% - Key parameter is d.TMAX: try with large values and then reduce if fails to converge

% rev: October 30, 2015
% rev: December, 2017
% rev: March, 2018 (Simone Ferrari)

clear all, close all, clc

SLASH='\'; 

% 1) Load magnetic model
load LastPath
[filename, pathname] = uigetfile([pathname 'fdfq*.mat'],'Load fdfq_idiq file');
filenamewithpath = [pathname filename];
load(filenamewithpath);
save LastPath pathname
% Id_ModMgn = Id; Iq_ModMgn = Iq;
% Fd_ModMgn = Fd; Fq_ModMgn = Fq;


% 2) Setup

prompt={'Motor name', ...
        'Pole pair number', ...
        'Axes type (SR or PM)', ...
        'Number of three phase sets',...
        'Max line-to-line peak voltage [V]', ...
        'Max phase peak current [A]',...
        'Min speed [rpm]', ...
        'Max speed [rpm]', ...
        'Min torque in the plot [Nm]',...
        'Max torque in the plot [Nm]',...
        'Number of speed level',...
        'Number of torque level',...
        'Reference temperature [°C]',...
        'Output temperature [°C]',...
        'Stator phase resistance [Ohm]',...
        'Iron Loss Model? (1=yes/0=no)',...
        'Skin Effect Model? (1=yes/0=no)',...
        'Mechanical Loss Model?'};
name='Input';
numlines=1;
defaultanswer={ 'MotName',...
                '2',...
                'SR',...
                '1',...
                '566',...
                '8.5',...
                '0',...
                '3000',...
                '0',...
                '14',...
                '31',...
                '15',...
                '120',...
                '120',...
                '4.6',...
                '0',...
                '0',...
                '[0]'};
setup=inputdlg(prompt,name,numlines,defaultanswer);

d.motorName = setup{1};
d.p         = eval(setup{2});
d.motorType = setup{3};
d.n3phase   = eval(setup{4});
d.Vbus      = eval(setup{5});
d.Imax      = eval(setup{6});
d.velmin    = eval(setup{7});
d.velmax    = eval(setup{8});
d.Tmin      = eval(setup{9});
d.Tmax      = eval(setup{10});
d.ns        = eval(setup{11});
d.nt        = eval(setup{12});
d.temp0     = eval(setup{13});
d.temp      = eval(setup{14});
d.Rs        = eval(setup{15});

flagIron    = eval(setup{16});
flagSkin    = eval(setup{17});
d.pMechLoss = eval(setup{18});

% 3) Load loss model
if flagIron
    [filename_loss, pathname_loss] = uigetfile([pathname SLASH '.mat'],'Load Iron Loss Model');
    d.IronLossModel=loadIronLossModel([pathname_loss filename_loss]);
else
    d.IronLossModel=loadIronLossModel('0');
end

% 4) Skin effect factor for stator Joule loss
if flagSkin
    [filename_slot, pathname_slot] = uigetfile([pathname SLASH '.mat'],'Load results from slot skin effect analysis');
    d.SkinEffModel=loadSkinEffectModel([pathname_slot filename_slot]);
else
    d.SkinEffModel=loadSkinEffectModel('0');
end

%d.velbase = d.velmax;

DatiOpt=PuntiDiOttimo(d,filenamewithpath);

DatiOpt.temperature=d.temp;


% print figures
StampaMaxTW;

ButtonName = questdlg('SAVE RESULTS?', 'Figure Export', 'No');
if strcmp(ButtonName,'Yes')
    % Save the results
    %OutputFolder=[pathname 'Output' datestr(now,30)];
    %OutputFolder=[pathname 'TwMap_' int2str(d.Tmax) 'Nm_' int2str(d.temp) 'C'];
    OutputFolder=[pathname 'TwMap_' d.motorName];
    mkdir(OutputFolder)
    if isoctave()  
        filenamewithpath=strcat(OutputFolder, SLASH, 'TwMap_', int2str(d.Tmax), 'Nm.mat');
        save ('-mat-binary', filenamewithpath);
        clear filenamewithpath
    else
        filenamewithpath=[OutputFolder SLASH 'TwMap_' int2str(d.Tmax) 'Nm.mat'];
        save(filenamewithpath,'DatiOpt','d')
    end
    
    % salva le figure
    %     [SUCCESS,MESSAGE,MESSAGEID] = mkdir([OutputFolder SLASH 'Figure_MaxTW']);
    for ii = 1:get(gcf,'Number') %H_last
        if isoctave()
            fig_name=strcat(OutputFolder, SLASH, 'fig_', num2str(ii));
            hgsave(ii,[fig_name]);
            clear fig_name
        else
            saveas(ii,[OutputFolder SLASH 'fig_' num2str(ii)],'fig');
        end        
    end

%     ButtonName = questdlg('SAVE RESULTS?', 'Excel export', 'No');
%     if strcmp(ButtonName,'Yes')
%         export_mapTw_excel;
%     end
    
    
end




