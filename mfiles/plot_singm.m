% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.


% interp the flux linkage maps over a very dense grid
id = F_map.Id; iq = F_map.Iq;
Fd0 = F_map.Fd; Fq0 = F_map.Fq;
T = F_map.T;
dT = F_map.dT;

i_d=linspace(id(1),id(end),n2);
i_q=linspace(iq(1),iq(end),n2);

[Id,Iq]=meshgrid(i_d,i_q);

% refine maps to the 256 x 256 standard resolution
Fd = interp2(idvect,iqvect,Fd0,Id,Iq,'cubic');
Fq = interp2(idvect,iqvect,Fq0,Id,Iq,'cubic');

T = interp2(idvect,iqvect,T,Id,Iq,'cubic');
T = T*klength;

dT = interp2(idvect,iqvect,dT,Id,Iq,'cubic');
dT = dT*klength;

%% rewind
Id=Id/kturns;
Iq=Iq/kturns;
Fd=Fd*kturns;
Fq=Fq*kturns;

%% adapt the stack lenght
Fd=Fd*klength;
Fq=Fq*klength;

% %% add end-connections term
% Fd = Fd + Lld * Id;
% Fq = Fq + Llq * Iq;

save ([NewDir 'fdfq_idiq_n' num2str(n2) '.mat'],'Fd','Fq','Id','Iq');
save ([NewDir 'fdfq_idiq_n' num2str(n2) '.mat'],'T','-append');
save ([NewDir 'fdfq_idiq_n' num2str(n2) '.mat'],'dT','-append');

% flux maps
figure
plot(Id(1,:),Fd([1 end],:),F_map.Id(1,:),F_map.Fd([1 end],:),'kx'), grid on, hold on
plot(Iq(:,1),Fq(:,[1 end]),F_map.Iq(:,1),F_map.Fq(:,[1 end]),'kx'), 
xlabel('id,iq [A]'), ylabel('\lambda_d, \lambda_q [Vs]'), %zlabel('\lambda_d')
saveas(gcf,[NewDir 'Curves_' strrep(filemot,'.mat','.fig')])

figure
surfc(Id,Iq,Fd), grid on, xlabel('id'), ylabel('iq'), zlabel('\lambda_d')
if not(kturns == 1)
    title(['Rewind factor = ' num2str(kturns)])
end
saveas(gcf,[NewDir 'Fdsurf_' strrep(filemot,'.mat','.fig')])

figure
surfc(Id,Iq,Fq), grid on, xlabel('id'), ylabel('iq'), zlabel('\lambda_q')
if not(kturns == 1)
    title(['Rewind factor = ' num2str(kturns)])
end
saveas(gcf,[NewDir 'Fqsurf_' strrep(filemot,'.mat','.fig')])


%% TORQUE MAP
figure
surf(Id,Iq,abs(T)), grid on, xlabel('id [A]'), ylabel('iq [A]'), zlabel('Torque [Nm]')
if not(kturns == 1)
    title(['Rewind factor = ' num2str(kturns)])
end
saveas(gcf,[NewDir 'Tsurf_' strrep(filemot,'.mat','.fig')])

figure
surf(Id,Iq,dT), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('Torque ripple [Nm]')
if not(kturns == 1)
    title(['Rewind factor = ' num2str(kturns)])
end
saveas(gcf,[NewDir 'dTdsurf_' strrep(filemot,'.mat','.fig')])

