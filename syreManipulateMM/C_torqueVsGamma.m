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

%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
% C_torqueVsGamma.m
% evaluates torque versus current phase angle gamma

% input:    fdfq_idiq_n256 (direct model) and ReadParameters.m
% output:   figures, saved in pathname

% rev: July 22, 2015

clear all, close all, addpath mfiles %addpath C:\Matlab_functions\
load LastPath
[FileName,pathname] = uigetfile([pathname '/*_n*.mat'],'LOAD DATA')
load([pathname FileName]); % run([pathname 'ReadParameters'])
save LastPath pathname

if not(exist('ID'))
    ID = Id;
    IQ = Iq;
end

if exist('T')
    TI = abs(T);
    p = round(abs(2/3 * T(2,2)./(Fd(2,2) .* Iq(2,2) - Fq(2,2) .* Id(2,2))));   % pole pairs reconstruction
else
    run([pathname 'ReadParameters']);
    TI = 3/2 * p * (Fd .* Iq - Fq .* Id);   % torque
end

Imax_interp = max(max(max(abs(Id))),max(max(abs(Iq))));
ILevel = Imax_interp * [0.25 0.5 0.75 1.00];    % current amplitude levels

gamma = linspace(0,90,20);

id = ILevel' * cosd(gamma);
iq = ILevel' * sind(gamma);

if sum(Id(:,1)) < 0
        axes_type = 'PM';  % PM style
    else
        axes_type = 'SR';  % SyR style
end

if axes_type == 'PM'
    id = ILevel' * -sind(gamma);
    iq = ILevel' * cosd(gamma);
end
    
t = interp2(ID,IQ,TI,id,iq);
fd = interp2(ID,IQ,Fd,id,iq);
fq = interp2(ID,IQ,Fq,id,iq);
PF = abs(sin(atan2(iq,id) - atan2(fq,fd)));

%% Output

% result folder
resFolder='plotVSgamma\';
mkdir(pathname,resFolder);


% torque
figure()
figSetting
hl=plot(gamma,abs(t),'-x');
xlabel('$\gamma$ [$^\circ$]')
ylabel('$T$ [Nm]')
for ii=1:length(ILevel)
    set(hl(ii),'DisplayName',[num2str(ILevel(ii),3) ' A'])
end
legend('show','Location','NorthWest')
set(gca,'XLim',[0 90],'XTick',0:15:90)
saveas(gcf,[pathname resFolder 'torqueVSgamma.fig']);

% power factor
figure()
figSetting
hl=plot(gamma,PF,'-x');
xlabel('$\gamma$ [$^\circ$]')
ylabel('$cos \varphi$')
for ii=1:length(ILevel)
    set(hl(ii),'DisplayName',[num2str(ILevel(ii),3) ' A'])
end
legend('show','Location','NorthWest')
set(gca,'XLim',[0 90],'XTick',0:15:90)
saveas(gcf,[pathname resFolder 'pfVSgamma.fig']);

% flux linkage
figure()
figSetting
hl=plot(gamma,abs(fd+j*fq),'-x');
xlabel('$\gamma$ [$^\circ$]')
ylabel('$\lambda$ [Vs]')
for ii=1:length(ILevel)
    set(hl(ii),'DisplayName',[num2str(ILevel(ii),3) ' A'])
end
legend('show','Location','SouthWest')
set(gca,'XLim',[0 90],'XTick',0:15:90)
saveas(gcf,[pathname resFolder 'fluxVSgamma.fig']);


% delta angle
figure()
figSetting
hl=plot(gamma,angle(fd+j*fq)*180/pi,'-x');
xlabel('$\gamma$ [$^\circ$]')
ylabel('$\delta$ [$^\circ$]')
for ii=1:length(ILevel)
    set(hl(ii),'DisplayName',[num2str(ILevel(ii),3) ' A'])
end
legend('show','Location','NorthWest')
set(gca,'XLim',[0 90],'XTick',0:15:90)
%set(gca,'YLim',[0 90],'YTick',0:15:90)
saveas(gcf,[pathname resFolder 'deltaVSgamma.fig']);

% gamma1 = 55;
% f1 = interp1(gamma,abs(fd(1,:)+j*fq(1,:)),gamma1)
% t1 = interp1(gamma,abs(t(1,:)),gamma1)
% n1 = 2800
% n2 = 8000
% gamma2 = interp1(abs(fd(1,2:end-1)+j*fq(1,2:end-1)),gamma(2:end-1),f1*n1/n2)
% 


