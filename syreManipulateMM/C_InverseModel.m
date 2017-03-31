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

%% C_InverseModel.m

% GP - 01 02 2015
% input:    direct model fdfq_idiq_n256.mat 
% output:   inverse model ID_dato_fd_fq, IQ_dato_fd_fq

% output files are in new subfolder AOA 

close all, clear all; load LastPath; %addpath m

% load motor data
[FILENAME, pathname, FILTERINDEX] = uigetfile([pathname '/*_n*.mat'], 'LOAD DATA');
save('ultimo','pathname');
load([pathname FILENAME]);
prompt={'p: # of pole pairs?'};
name='prompt for p';
numlines=1;
defaultanswer={'2'};
answer=inputdlg(prompt,name,numlines,defaultanswer);
p = eval(answer{1});
[success,message,messageid] = mkdir(pathname,'AOA')
pathname = [pathname,'AOA\'];

C_IdIq_FqFd;

figure
temp = max(max(abs(FF)));
[c,h]=contour(FD,FQ,FF,temp * [0.1:0.1:0.9]); clabel(c,h,'Color','k'); hold on
temp = max(max(abs(TF)));
[c,h]=contour(FD,FQ,TF); clabel(c,h);
% [c,h]=contour(FD,FQ,TF,temp * [-0.9:0.1:-0.1 -0.1:0.01:0.1 0.1:0.1:0.9]); clabel(c,h);
% [c,h]=contour(FD,FQ,IF,'k'); clabel(c,h);
axis equal, grid on
xlabel('fd [Vs]'), ylabel('fq [Vs]'), title('flux linkage and torque contours')

if exist('fd_KtMax')
    plot(fd_KtMax,fq_KtMax,'k','LineWidth',2)
    plot(fd_KvMax,fq_KvMax,'k','LineWidth',2)
    saveas(gcf,[pathname 'MTPA-MTPV - flux']);
end