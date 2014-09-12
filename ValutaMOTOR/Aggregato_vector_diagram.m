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


clear all, close all
Foldername='ValutaMOTOR';
PATHNAME=cd;
pathname=PATHNAME(1:end-length(Foldername));
addpath([pathname,'Matlab_Functions'])
load ultimo.mat
% save('ultimo','pathname');
if pathname == 0
    pathname = [cd '\'];
end
num_casi=3;

kl=1; kspire=1;
color={'r','g','b'};
legenda={};
for i=1:num_casi
[filename pathname] = uigetfile([pathname '*_*_*.mat'],'pick a file');
load([pathname filename(1:end-4)]);
save('ultimo','pathname');

t = out.SOL(:,6)*kl;

id= mean(out.SOL(:,2))/kspire;
iq= mean(out.SOL(:,3))/kspire;

fd = mean(out.SOL(:,4))*kl*kspire;
fq = mean(out.SOL(:,5))*kl*kspire;

I=abs(id+1j*iq);
gamma=atan(iq/id)*180/pi;

F=abs(fd+1j*fq);
delta=atan(iq/id)*180/pi;
 filename(filename=='_')='-';
 legenda{i}=filename(1:end-4);
figure(1); plot([0,id]',[0,iq]',color{i},'LineWidth',2);hold on;
figure(2); plot([0,fd]',[0,fq]',color{i},'LineWidth',2);hold on;

end
 hold off;
figure(1);axis square; legend(legenda);grid on; xlabel('id'); ylabel('iq');
figure(2);axis square; legend(legenda);grid on; xlabel('fd'); ylabel('fq');