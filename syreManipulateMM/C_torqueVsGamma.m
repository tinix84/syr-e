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
    id = ILevel * -sind(gamma);
    iq = ILevel * cosd(gamma);
end
    
t = interp2(ID,IQ,TI,id,iq);
fd = interp2(ID,IQ,Fd,id,iq);
fq = interp2(ID,IQ,Fq,id,iq);
PF = sind(atand(iq./id) - atand(fq./fd));

% torque
figure
plot(gamma,abs(t),'-x'), grid on
xlabel('current phase angle [deg]')
ylabel('Torque [Nm]')
legend([num2str(ILevel,3) ' A'])
%adapt_figure(0.4), adapt_figure_fonts('Times New Roman',12,10)
xlim([0 90]), saveas(gcf,[pathname 'coppia_su_gamma'])

% power factor
figure
plot(gamma,PF,'-x'), grid on
xlabel('current phase angle [deg]')
ylabel('Power Factor')
legend([num2str(ILevel,3) ' A'])
%adapt_figure(0.4), adapt_figure_fonts('Times New Roman',12,10)
xlim([0 90]), saveas(gcf,[pathname 'PF_su_gamma'])

% flux linkage
figure
plot(gamma,abs(fd+j*fq),'-x'), grid on
xlabel('current phase angle [deg]')
ylabel('Flux Linkage [Vs]')
legend([num2str(ILevel,3)])
%adapt_figure(0.4), adapt_figure_fonts('Times New Roman',12,10)
xlim([0 90]), saveas(gcf,[pathname 'flux_su_gamma'])

% delta angle
figure
plot(gamma,angle(fd+j*fq)*180/pi,'-x'), grid on
xlabel('current phase angle [deg]')
ylabel('\delta [deg]')
legend([num2str(ILevel,3)])
%adapt_figure(0.4), adapt_figure_fonts('Times New Roman',12,10)
xlim([0 90]), saveas(gcf,[pathname 'coppia_su_delta'])

% gamma1 = 55;
% f1 = interp1(gamma,abs(fd(1,:)+j*fq(1,:)),gamma1)
% t1 = interp1(gamma,abs(t(1,:)),gamma1)
% n1 = 2800
% n2 = 8000
% gamma2 = interp1(abs(fd(1,2:end-1)+j*fq(1,2:end-1)),gamma(2:end-1),f1*n1/n2)
% 


