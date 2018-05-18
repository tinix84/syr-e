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

function C_MtpaMtpvLut(pathname,filename)
% evaluates MTPA and MTPV in id,iq and in fd,fq

% input:    1) fdfq_idiq_n256 (direct model)
%           2) ReadParameters.m
% output:   1) MTPA (called kt max), MTPV (called kv max)
%           2) tables printed in tables_MTPA.txt

% output files are in new subfolder AOA, created into the folder where
% machine data resides

% Definitions:
% direct model: Id(n,n),Iq(n,n) -> Fd(n,n), Fq(n,n)
% inverse model (evaluated here): fd(1,m),fq(1,n) -> ID_dato_fd_fq (n,m), IQ_dato_fd_fq (n,m)

% release: July 22, 2015


%% old input
% close all, clear all; addpath mfiles, load LastPath
%
% % load the two motor data: mat file and ReadParameters.m
% % 1) flux linkage map
% [FILENAME, pathname, FILTERINDEX] = uigetfile([pathname '/*_n*.mat'], 'LOAD DATA');
% % 2) associated ReadParameters.m
% load([pathname FILENAME]); %run([pathname 'ReadParameters']);
% save LastPath pathname
%
% [success,message,messageid] = mkdir(pathname,'AOA')
% GoAhead = 'Yes';
% if not(isempty(message))
%     GoAhead = questdlg('Existing folder: proceed anyway?', ...
%         'WARNING!', ...
%         'Yes', 'No', 'No');
% end

%% new input

close all

if nargin()<2
    load LastPath
    [filename, pathname, ~] = uigetfile([pathname '/*_n*.mat'], 'LOAD DATA');
    save LastPath pathname
    
    [success,message,messageid] = mkdir(pathname,'AOA');
    GoAhead = 'Yes';
    %     if not(isempty(message))
    %         GoAhead = questdlg('Existing folder: proceed anyway?', ...
    %             'WARNING!', ...
    %             'Yes', 'No', 'No');
    %     end
else
    [success,message,messageid] = mkdir(pathname,'AOA');
    if ~isempty(message)
        warning('Existing folder: the old data will be lost');
    end
end

load([pathname filename]);

%plot_defaults;

% if isempty(axes_type)
if sum(Id(:,1)) < 0
    axes_type = 'PM';  % SPM style
else
    axes_type = 'SR';  % IPM style
end

% end

if(0)
    % change mesh resolution
    npoints = 32;
    i_d=linspace(Id(1),Id(end),npoints); i_q=linspace(Iq(1),Iq(end),npoints);
    Id0 = Id; Iq0 = Iq; Ld0 = Fd; Lq0 = Fq;
    [Id,Iq]=meshgrid(i_d,i_q);
    Fd = interp2(Id0,Iq0,Ld0,Id,Iq,'cubic');
    Fq = interp2(Id0,Iq0,Lq0,Id,Iq,'cubic');
end

if (0)
    % plot the initial magnetic model
    plot_maps(pathname,filename);
end

if (0)
    % keep only the data of positive dq currents
    
    answ=inputdlg('axes type (SR or PM)?','Axes type',1,{axes_type});
    axes_type=answ{1};
    
    Id0=Id;
    Iq0=Iq;
    Fd0=Fd;
    Fq0=Fq;
    T0=T;
    dTpp0=dTpp;
    dT0=dT;
    
    Id=linspace(0,max(max(Id)),256);
    Iq=linspace(0,max(max(Iq)),256);
    [Id,Iq]=meshgrid(Id,Iq);
    Fd=interp2(Id0,Iq0,Fd0,Id,Iq);
    Fq=interp2(Id0,Iq0,Fq0,Id,Iq);
    T=interp2(Id0,Iq0,T0,Id,Iq);
    dT=interp2(Id0,Iq0,dT0,Id,Iq);
    dTpp=interp2(Id0,Iq0,dTpp0,Id,Iq);
end


% Performance maps
id = Id(1,:); iq = Iq(:,1)';
if exist('T')
    TI = abs(T);
    if ~exist('n3phase')
        n3phase=1;
        warning('Number of 3 phase sets equal to one as default setting')
    end
    p = round(abs(2/3 * T(2,2)./(Fd(2,2) .* Iq(2,2) - Fq(2,2) .* Id(2,2))/n3phase));   %AS pole pairs reconstruction
else
    if exist([pathname 'ReadParameters'],'file')
        run([pathname 'ReadParameters']);
    else
        prompt={'p: # of pole pairs?'};
        name='prompt for p';
        numlines=1;
        defaultanswer={'2'};
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        p = eval(answer{1});
    end
    TI = 3/2 * p * (Fd .* Iq - Fq .* Id);   % torque
end
FI = sqrt(Fd.^2 + Fq.^2);               % flux linkage amplitude
II = sqrt(Id.^2 + Iq.^2);               % current amplitude
IPF = cos(atan(Fq./Fd)-atan(-Id./Iq));  % internal PF

if axes_type == 'PM'
    fm = interp2(Id,Iq,FI,-eps,eps);    % flux linkage at open circuit
else
    fm = interp2(Id,Iq,FI,eps,eps);
end

if (0)
    % Inductances (no file saved)
    C_ldlq_idiq
end

pathname = [pathname 'AOA\'];
if (0)
    % Inverse Magnetic Model, saves IdIq_FdFq.mat
    C_IdIq_FqFd
    
    [FD,FQ]=meshgrid(fd,fq);
    T = 3/2 * p * n3phase * (FD .* IQ_dato_fd_fq - FQ .* ID_dato_fd_fq); %AS
    F = sqrt(FD.^2 + FQ.^2);
    delta = zeros(size(FD));
    delta(FD ~= 0) = atan(FQ(FD ~= 0)./FD(FD ~= 0));
    I = sqrt(ID_dato_fd_fq.^2 + IQ_dato_fd_fq.^2);
    ID_dato_fd_fq(ID_dato_fd_fq == 0) = eps;
    argI = atan(IQ_dato_fd_fq./ID_dato_fd_fq);
    fi = pi/2 + delta - argI;    % PF angle
    PF = cos(fi);
    
    
    % flux coordinates - contours
    CurveIsoT=contourc(fd,fq,T,TLevel);
    CurveIsoI=contourc(fd,fq,I,ILevel);
    CurveIsoPF=contourc(fd,fq,PF,PFLevel);
    
    % MTPV
    KvMax_Locus_FdFq;   % in id,iq
    
end
pathname = strrep(pathname,'AOA\','');

if (0)
    plot_all_surfaces;
end

% CONTOURS, used for MTPA and MTPV evaluation
TLevel = min(max(max(TI)),max(max(TI)))*(0.01:0.01:1);    % torque levels
Imax_interp = max(max(abs(id)),max(abs(iq)));
ILevel = Imax_interp *(0.01:0.01:1);    % current amplitude levels
PFLevel = 0.5:0.1:0.9;                  % PF levels

% current coordinates - contours
CurveIsoTI=contourc(id,iq,TI,TLevel);
CurveIsoII=contourc(id,iq,II,ILevel);
CurveIsoTpF = contourc(id,iq,TI./FI,200);

pathname = [pathname 'AOA\'];

% MTPA
KtMax_Locus_idiq;   % in id,iq
KtMax_idiq2fdfq;    % in fd,fq
F_KtMax = abs(fd_KtMax + 1i* fq_KtMax);

% MTPV
% KvMax_Locus_FdFq;   % in id,iq
KvMax_Locus_idiq;   % in fd,fq

pathname = strrep(pathname,'AOA\','');

% Performance Curves
figure
if ~isoctave()
   figSetting
end    
plot(I_KtMax,T_KtMax), grid on,
if exist('dTpp_KtMax')
    hold on
    plot(I_KtMax,[(T_KtMax+dTpp_KtMax/2);(T_KtMax-dTpp_KtMax/2)],'r');
end
set(gca,'PlotBoxAspectRatio',[1 0.6 1]);
xlabel('Phase current [A] pk'), ylabel('Torque [Nm]')
title('Torque vs pk current, MPTA')
h=gcf(); %AS
if isoctave()
    fig_name=strcat(pathname, 'AOA\torque_current');
    hgsave(h,[fig_name]);
    clear fig_name
else
    saveas(gcf,[pathname 'AOA\torque_current']);
end

figure
if ~isoctave()
   figSetting
end  
% gamma_KtMax = atand(-id_KtMax./iq_KtMax);
gamma_KtMax = atand(iq_KtMax./id_KtMax);
plot(I_KtMax,gamma_KtMax), grid on
set(gca,'PlotBoxAspectRatio',[1 0.6 1]);
xlabel('Phase current [A] pk'), ylabel('MTPA current angle [deg]')
title('$$\gamma_{MPTA}(I)$$')
h=gcf(); %AS
if isoctave()
    fig_name=strcat(pathname, 'AOA\gamma_current');
    hgsave(h,[fig_name]);
    clear fig_name
else
    saveas(gcf,[pathname 'AOA\gamma_current']);
end

figure
if ~isoctave()
   figSetting
end  
plot(T_KtMax,F_KtMax), grid on
set(gca,'PlotBoxAspectRatio',[1 0.6 1]);
xlabel('Torque [Nm]'), ylabel('Flux amplitude [Vs]')
title('Torque vs pk flux, MPTA')
h=gcf(); %AS
if isoctave()
    fig_name=strcat(pathname, 'AOA\flux_torque');
    hgsave(h,[fig_name]);
    clear fig_name
else
    saveas(gcf,[pathname 'AOA\flux_torque']);
end

cut_point = 7;
figure
if ~isoctave()
   figSetting
end  
plot(I_KtMax(cut_point:end),T_KtMax(cut_point:end)./I_KtMax(cut_point:end)), grid on
set(gca,'PlotBoxAspectRatio',[1 0.6 1]);
xlabel('Phase current [A] pk'), ylabel('kt [Nm/Apk]')
title('kt vs pk current, MPTA')
h=gcf(); %AS
if isoctave()
    fig_name=strcat(pathname, 'AOA\kt_current');
    hgsave(h,[fig_name]);
    clear fig_name
else
    saveas(gcf,[pathname 'AOA\kt_current']);
end


% contours
figure
if ~isoctave()
   figSetting
end  
[c,h]=contour(id,iq,TI,'k'); clabel(c,h,'LabelSpacing',200), hold on
[c,h]=contour(id,iq,II,'k'); clabel(c,h)
plot(id_KtMax(1:end),iq_KtMax(1:end),'-b','LineWidth',2),
plot(id_KvMax,iq_KvMax,'-b','LineWidth',2),
grid on, xlabel('$$i_d$$ [A]'),ylabel('$$i_q$$ [A]'), hold off
axis equal, axis([min(id) max(id) min(iq) max(iq)])
title('Control curves in id iq coordinates')
legend('[Nm]','[A]','MTPA','MTPV')
% adapt_figure_fonts('Times New Roman',18,14)
h=gcf(); %AS
if isoctave()
    fig_name=strcat(pathname, 'AOA\MTPAMTPV(id,iq)');
    hgsave(h,[fig_name]);
    clear fig_name
else
    saveas(gcf,[pathname 'AOA\MTPAMTPV(id,iq)']);
end

if exist('F')
    figure
    if ~isoctave()
       figSetting
    end  
    temp = max(max(abs(F)));
    [c,h]=contour(FD,FQ,F,temp * [0.1:0.1:0.9]);
    clabel(c,h,'Color','k'); hold on
    temp = max(max(abs(TF)));
    [c,h]=contour(FD,FQ,TF,temp * (-0.9:0.05:0.9)); clabel(c,h);
    [c,h]=contour(FD,FQ,IF,'k'); clabel(c,h);
    axis equal, grid on
    xlabel('fd [Vs]'),ylabel('fq [Vs]')
    title('AOA in fd,fq coordinates')
    plot(fd_KtMax,fq_KtMax,'-b','LineWidth',2)   % MTPA
    plot(fd_KvMax,fq_KvMax,'-b','LineWidth',2)   % MTPV
    h=gcf(); %AS
    if isoctave()
        fig_name=strcat(pathname, 'AOA\MTPAMTPV(fd,fq)');
        hgsave(h,[fig_name]);
        clear fig_name
    else
        saveas(gcf,[pathname 'AOA\MTPAMTPV(fd,fq)']);
    end
end

% print MTPA tables into a txt file
T_KtMax(1) = 0;

% LUT
m = 1;  % # of lines
n = 20; % # of columns (table size is 1 x n)

Tmax = T_KtMax(end);
step = Tmax/n;
T_set = 0:step:Tmax;

id_set = interp1(T_KtMax,id_KtMax,T_set); id_set(1) = 0;
iq_set = interp1(T_KtMax,iq_KtMax,T_set); iq_set(1) = 0;
id_set(1) = 0;
iq_set(1) = 0;
fd_set = interp1(T_KtMax,fd_KtMax,T_set);
fq_set = interp1(T_KtMax,fq_KtMax,T_set);
f_set = interp1(T_KtMax,abs(fd_KtMax + j*fq_KtMax),T_set);


% print txt file (MTPA)
fid = fopen([pathname 'tablesMTPA.txt'],'w');
% fprintf(fid,'//SIGLA MOTORE: %s\n',motor_name);
fprintf(fid,['//' date '\n']);
fprintf(fid,'float TMIN    = 0;\n');
fprintf(fid,'float TMAX    = %4.3f; //Nm\n',Tmax);
fprintf(fid,'float DT      = %4.4f; //Nm\n',step);
fprintf(fid,'float INV_DT  = %4.4f; //Nm^-1\n',1/step);

StampaVarg(fid,id_set,m,n+1,'ID_REF','//MTPA - id','%6.3f')
StampaVarg(fid,iq_set,m,n+1,'IQ_REF','//MTPA - iq','%6.3f')
StampaVarg(fid,fd_set,m,n+1,'FD_REF','//MTPA - fd','%6.3f')
StampaVarg(fid,fq_set,m,n+1,'FQ_REF','//MTPA - fq','%6.3f')
StampaVarg(fid,f_set,m,n+1,'F_REF','//MTPA - flux amplitude','%6.3f')

fclose(fid);

% debug
figure
if ~isoctave()
   figSetting
end  
set(gca,'PlotBoxAspectRatio',[1 0.6 1]);
plot(T_KtMax,id_KtMax), hold on, grid on,
plot(T_KtMax,iq_KtMax,'r'),
plot(T_set,id_set,'xk'), plot(T_set,iq_set,'xk'), hold off,
set(gca,'PlotBoxAspectRatio',[1 0.6 1]);
legend('id setpoint','iq setpoint','LUT','LUT');
xlabel('$$T\_set [Nm]$$'), ylabel('[A] pk')
h=gcf(); %AS
if isoctave()
    fig_name=strcat(pathname, 'AOA\tablesMTPA');
    hgsave(h,[fig_name]);
    clear fig_name
else
    saveas(gcf,[pathname 'AOA\tablesMTPA']);
end

edit([pathname '\tablesMTPA.txt'])