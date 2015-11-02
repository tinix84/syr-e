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

% ...

close all, clear all; addpath mfiles, load LastPath

% load the two motor data: mat file and ReadParameters.m
% 1) flux linkage map
[FILENAME, pathname, FILTERINDEX] = uigetfile([pathname '/*_n*.mat'], 'LOAD DATA');
% 2) associated ReadParameters.m
load([pathname FILENAME]); %run([pathname 'ReadParameters']);
save LastPath pathname

[success,message,messageid] = mkdir(pathname,'AOA')
GoAhead = 'Yes';
if not(isempty(message))
    GoAhead = questdlg('Existing folder: proceed anyway?', ...
        'WARNING!', ...
        'Yes', 'No', 'No');
end

% if isempty(axes_type)
    if sum(Id(:,1)) < 0
        axes_type = 'PM';  % SPM style
    else
        axes_type = 'SR';  % IPM style
    end
% end

% plot the initial magnetic model
plot_maps(pathname,FILENAME);
 
% LUTs for flux observer
id_tab_min = min(min(Id));
id_tab_max = max(max(Id));
iq_tab_min = min(min(Iq));
iq_tab_max = max(max(Iq));

% LUT dimension
m = 9;    % rows
n = 21;   % columns

% Fd table: current steps
Didd = (id_tab_max-id_tab_min)/(n-1);
Diqd = (iq_tab_max-iq_tab_min)/(m-1);
% Fq table: current steps
Diqq = (iq_tab_max-iq_tab_min)/(n-1);
Didq = (id_tab_max-id_tab_min)/(m-1);

[idd,iqd]=meshgrid(linspace(id_tab_min,id_tab_max,n),linspace(iq_tab_min,iq_tab_max,m));
[idq,iqq]=meshgrid(linspace(id_tab_min,id_tab_max,m),linspace(iq_tab_min,iq_tab_max,n));

fd=interp2(Id,Iq,Fd,idd,iqd);
fq=interp2(Id,Iq,Fq,idq,iqq);
fq=fq';

% print to FluxTables.txt
fid = fopen([pathname 'FluxTables.txt'],'w');
% fprintf(fid,'//MOTOR NAME: %s\n',motor_name);
fprintf(fid,['//' date '\n']);
fprintf(fid,['float  ID_TAB_MIN = ' num2str(id_tab_min) ' ;\r\n']);
fprintf(fid,['float  IQ_TAB_MIN = ' num2str(iq_tab_min) ' ;\r\n']);
fprintf(fid,['float  ID_TAB_MAX = ' num2str(id_tab_max) ' ;\r\n']);
fprintf(fid,['float  IQ_TAB_MAX = ' num2str(iq_tab_max) ' ;\r\n']);
fprintf(fid,['float  DIDD = ' num2str(Didd,4) ' ;\r\n']);
fprintf(fid,['float  DIQD = ' num2str(Diqd,4) ' ;\r\n']);
fprintf(fid,['float  DIQQ = ' num2str(Diqq,4) ' ;\r\n']);
fprintf(fid,['float  DIDQ = ' num2str(Didq,4) ' ;\r\n']);
fprintf(fid,['float  INV_DIDD = ' num2str(1/Didd,4) ' ;\r\n']);
fprintf(fid,['float  INV_DIQD = ' num2str(1/Diqd,4) ' ;\r\n']);
fprintf(fid,['float  INV_DIQQ = ' num2str(1/Diqq,4) ' ;\r\n']);
fprintf(fid,['float  INV_DIDQ = ' num2str(1/Didq,4) ' ;\r\n']);
% fd(:,1)=0*fd(:,1);  % Fd @ id=0
StampaVarg(fid,fd',m,n,'FD_LUT','//Fluxd(iq,id)','%6.4f')
StampaVarg(fid,fq',m,n,'FQ_LUT','//Fluxq(id,iq)','%6.4f')
fclose(fid);

figure,
plot(idd',fd'), grid on, hold on
plot(iqq ,fq'),
plot(idd',fd','kx'), 
plot(iqq ,fq','kx'), hold off
xlabel('i_d, i_q [A]'), ylabel('\lambda_d \lambda_q [Vs]'), %title(motor_name)

saveas(gcf,[pathname 'FluxTables'])
edit([pathname 'FluxTables.txt'])