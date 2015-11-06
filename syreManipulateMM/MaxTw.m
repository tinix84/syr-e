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

clear all, close all, clc, addpath C:\Matlab_Functions, addpath mfiles
SLASH='\'; 

% 1) Load magnetic model
load LastPath
[filename, pathname] = uigetfile([pathname SLASH 'fdfq*.mat'],'IdIqLdLq');
filenamewithpath = [pathname filename];
load(filenamewithpath);
save LastPath pathname
Id_ModMgn = Id; Iq_ModMgn = Iq;
Fd_ModMgn = Fd; Fq_ModMgn = Fq;

% 2) Load loss model (can be same file)
[filename_loss, pathname] = uigetfile([pathname SLASH 'fdfq*.mat'],'LOSS MAP');
if ischar(filename_loss)
    load([pathname filename_loss]);
end
Id_Loss = Id; Iq_Loss = Iq;
Fd_Loss = Fd; Fq_Loss = Fq;

% 3) Load other inputs
filenamewithpath=[pathname 'ReadParameters.m'];
run(filenamewithpath)

% current region, common to both models
id_min = min([min(min(Id_ModMgn)),min(min(Id_Loss))]); id_max = max([max(max(Id_ModMgn)),max(max(Id_Loss))]);
iq_min = min([min(min(Iq_ModMgn)),min(min(Iq_Loss))]); iq_max = max([max(max(Iq_ModMgn)),max(max(Iq_Loss))]);

id_vett=linspace(id_min,id_max,256); iq_vett=linspace(iq_min,iq_max,256);
[Id,Iq]=meshgrid(id_vett,iq_vett);

% Flux linkages
Fd=interp2(Id_ModMgn,Iq_ModMgn,Fd_ModMgn,Id,Iq,'spline'); Fq=interp2(Id_ModMgn,Iq_ModMgn,Fq_ModMgn,Id,Iq,'spline');

% Loss
if exist('Pfes_c','var')
    Pfes_c=interp2(Id_Loss,Iq_Loss,Pfes_c,Id,Iq,'spline'); Pfes_h=interp2(Id_Loss,Iq_Loss,Pfes_h,Id,Iq,'spline');
    Pfer_c=interp2(Id_Loss,Iq_Loss,Pfer_c,Id,Iq,'spline'); Pfer_h=interp2(Id_Loss,Iq_Loss,Pfer_h,Id,Iq,'spline');
end
if exist('P_BarRot','var')
    P_BarRot=interp2(Id_Loss,Iq_Loss,P_BarRot,Id,Iq,'spline');
end
if exist('Ppm','var')
    Ppm=interp2(Id_Loss,Iq_Loss,Ppm,Id,Iq,'spline');
end

% add zeros into empty spots
setZero = ((round(Id)>max(max(Id_Loss)))+(Id<min(min(Id_Loss)))+(Iq>max(max(Iq_Loss)))+(Iq<min(min(Iq_Loss))))>0;
if exist('Pfes_c','var')
    Pfes_c(setZero)=0; Pfes_h(setZero)=0;
    Pfer_c(setZero)=0; Pfer_h(setZero)=0;
end
if exist('P_BarRot','var')
    P_BarRot(setZero)=0;
end
if exist('Ppm','var')
    Ppm(setZero)=0;
end

% Prompt for user input
prompt={'change # turns (N''/N)', ...
    'change stack length (L''/L)', ...
    'end conn. length p.u.', ...
    '# torque levels in the grid', ...
    '# speed levs in the grid', ...
    'motor type (PM or SR)','PM loss: p.u. of total due to axial segmentation'};
name='Input';
numlines=1;
defaultanswer={'1','1','0','16','20','PM','1'};
setup=inputdlg(prompt,name,numlines,defaultanswer);

% compute control trajectories
if exist('Pfes_c','var')
    % iron loss map is available
    d.velDim = velDim;
    DatiOpt = PuntiDiOttimo(d,setup,Id,Iq,Fd,Fq,Pfes_c,Pfer_c,Pfes_h,Pfer_h,Ppm,coeff_Pfe);
else
    % iron loss is extrapolated from single test point d.PfeDim
    DatiOpt = PuntiDiOttimo(d,setup,Id,Iq,Fd,Fq);
end

% print figures
StampaMaxTW;

ButtonName = questdlg('SAVE RESULTS?', 'title', 'No');
if strcmp(ButtonName,'Yes')
    % Save the results
    OutputFolder=[pathname 'Output' datestr(now,30)];
    mkdir(OutputFolder)
    filenamewithpath=[OutputFolder SLASH 'Output' datestr(now,30) '.mat'];
    save(filenamewithpath,'DatiOpt','d')
    
    % salva le figure
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir([OutputFolder SLASH 'Figure_MaxTW']);
    for ii = 1:H_last
        saveas(ii,[OutputFolder SLASH 'fig_' num2str(ii)],'fig');
    end
end

% print txt
fid = fopen([pathname 'tablesMTPAFluxWeak.txt'],'w');
% fprintf(fid,'//MOTOR NAME: %s\n',motor_name);
fprintf(fid,['//' date '\n']);
fprintf(fid,'float TMIN    = 0;\n');
fprintf(fid,'float TMAX    = %4.3f; //Nm\n',DatiOpt.Tmap(end));
fprintf(fid,'float DT      = %4.4f; //Nm\n',diff(DatiOpt.Tmap(1:2)));
fprintf(fid,'float INV_DT  = %4.4f; //Nm^-1\n',1/diff(DatiOpt.Tmap(1:2)));

fprintf(fid,'float NMIN    = 0;\n');
fprintf(fid,'float NMAX    = %5.1f; //Nm\n',DatiOpt.velmec(end));
fprintf(fid,'float DN      = %5.1f; //Nm\n',diff(DatiOpt.velmec(1:2)));
fprintf(fid,'float INV_DN  = %5.1f; //Nm^-1\n',1/diff(DatiOpt.velmec(1:2)));


% fd(:,1)=0*fd(:,1);  % Fd @ id=0

[m,n] = size(DatiOpt.IdMin);

% eliminate NaN in txt print
ID_LUT = zeros(m,n); IQ_LUT = zeros(m,n);
for i = 1:m
    temp = DatiOpt.IdMin(i,:);
    temp1 = temp(isfinite(temp));
    temp(isnan(temp)) = temp1(end);
    ID_LUT(i,:) = temp;
    
    temp = DatiOpt.IqMin(i,:);
    temp1 = temp(isfinite(temp));
    temp(isnan(temp)) = temp1(end);
    IQ_LUT(i,:) = temp;   
end

FD_LUT = interp2(Id,Iq,Fd,ID_LUT,IQ_LUT);

StampaVarg(fid,ID_LUT',m,n,'ID_LUT_FW','//id(T,n)','%6.4f')
StampaVarg(fid,FD_LUT',m,n,'FD_LUT_FW','//fd(T,n)','%6.4f')

StampaVarg(fid,IQ_LUT',m,n,'IQ_LUT_FW','//iq(T,n)','%6.4f')

fclose(fid);

figure,
plot(DatiOpt.Tmap,ID_LUT), grid on, hold on
plot(DatiOpt.Tmap,ID_LUT,'kx'),
xlabel('[Nm]'), ylabel('i_d [A]'), %title(motor_name)

figure,
plot(DatiOpt.Tmap,IQ_LUT), grid on, hold on
plot(DatiOpt.Tmap,IQ_LUT,'kx'),
xlabel('[Nm]'), ylabel('i_q [A]'), %title(motor_name)

figure,
plot(DatiOpt.Tmap,FD_LUT), grid on, hold on
plot(DatiOpt.Tmap,FD_LUT,'kx'),
xlabel('[Nm]'), ylabel('\lambda_d [Vs]'), %title(motor_name)

% saveas(gcf,[pathname 'tablesMTPAFluxWeak'])
edit([pathname 'tablesMTPAFluxWeak.txt'])

