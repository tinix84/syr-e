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
Fd = F_map.Fd; Fq = F_map.Fq;

T = F_map.T;
if isfield(F_map,'dT')
    dT = F_map.dT;
end
if isfield(F_map,'dTpp')
    dTpp = F_map.dTpp;
end
if isfield(F_map,'Pfe')
    Pfe = F_map.Pfe;
end
if isfield(F_map,'Pfes_h')
    Pfes_h = F_map.Pfes_h;
    Pfes_c = F_map.Pfes_c;
    Pfer_h = F_map.Pfer_h;
    Pfer_c = F_map.Pfer_c;
end
if isfield(F_map,'Ppm')
    Ppm = F_map.Ppm;
end

i_d=linspace(id(1),id(end),n2);
i_q=linspace(iq(1),iq(end),n2);

[Id,Iq]=meshgrid(i_d,i_q);

% refine maps to the 256 x 256 standard resolution
% Fd = interp2(idvect,iqvect,Fd,Id,Iq,'cubic')*klength;
% Fq = interp2(idvect,iqvect,Fq,Id,Iq,'cubic')*klength;
% T = interp2(idvect,iqvect,T,Id,Iq,'cubic')*klength;

Fd = interp2(F_map.Id(1,:),F_map.Iq(:,1)',Fd,Id,Iq,'cubic')*klength;
Fq = interp2(F_map.Id(1,:),F_map.Iq(:,1)',Fq,Id,Iq,'cubic')*klength;
T = interp2(F_map.Id(1,:),F_map.Iq(:,1)',T,Id,Iq,'cubic')*klength;

if isfield(F_map,'dT')
    dT = interp2(F_map.Id(1,:),F_map.Iq(:,1)',dT,Id,Iq,'cubic');
    dT = dT*klength;
end
if isfield(F_map,'dTpp')
    dTpp = interp2(F_map.Id(1,:),F_map.Iq(:,1)',dTpp,Id,Iq,'cubic');
    dTpp = dTpp*klength;
end
if isfield(F_map,'Pfe')
    Pfe = interp2(F_map.Id(1,:),F_map.Iq(:,1)',Pfe,Id,Iq,'cubic');
end
if isfield(F_map,'Pfes_h')
    Pfes_h = interp2(F_map.Id(1,:),F_map.Iq(:,1)',Pfes_h,Id,Iq,'cubic')*klength;
    Pfes_c = interp2(F_map.Id(1,:),F_map.Iq(:,1)',Pfes_c,Id,Iq,'cubic')*klength;
    Pfer_h = interp2(F_map.Id(1,:),F_map.Iq(:,1)',Pfer_h,Id,Iq,'cubic')*klength;
    Pfer_c = interp2(F_map.Id(1,:),F_map.Iq(:,1)',Pfer_c,Id,Iq,'cubic')*klength;
end
if isfield(F_map,'Ppm')
    Ppm = interp2(F_map.Id(1,:),F_map.Iq(:,1)',Ppm,Id,Iq,'cubic')*klength;
end

%% rewind
Id=Id/kturns;
Iq=Iq/kturns;
Fd=Fd*kturns;
Fq=Fq*kturns;

% %% add end-connections term
% Fd = Fd + Lld * Id;
% Fq = Fq + Llq * Iq;

save ([NewDir 'fdfq_idiq_n' num2str(n2) '.mat'],'Fd','Fq','Id','Iq');
save ([NewDir 'fdfq_idiq_n' num2str(n2) '.mat'],'T','-append');
if isfield(F_map,'dT')
    save ([NewDir 'fdfq_idiq_n' num2str(n2) '.mat'],'dT','-append');
end
if isfield(F_map,'dTpp')
    save ([NewDir 'fdfq_idiq_n' num2str(n2) '.mat'],'dTpp','-append');
end
if isfield(F_map,'Pfe')
    save ([NewDir 'fdfq_idiq_n' num2str(n2) '.mat'],'Pfe','-append');
end
if isfield(F_map,'Pfes_h')
    save ([NewDir 'fdfq_idiq_n' num2str(n2) '.mat'],'Pfes_h','Pfes_c','Pfer_h','Pfer_c','-append');
end
if isfield(F_map,'Ppm')
    save ([NewDir 'fdfq_idiq_n' num2str(n2) '.mat'],'Ppm','-append');
end

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

if isfield(F_map,'dT')
    figure
    surf(Id,Iq,dT), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('Torque ripple (std) [Nm]')
    if not(kturns == 1)
        title(['Rewind factor = ' num2str(kturns)])
    end
    saveas(gcf,[NewDir 'dTsurf_' strrep(filemot,'.mat','.fig')])
end
if isfield(F_map,'dTpp')
    figure
    surf(Id,Iq,dTpp), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('Torque ripple (pk-pk) [Nm]')
    if not(kturns == 1)
        title(['Rewind factor = ' num2str(kturns)])
    end
    saveas(gcf,[NewDir 'dTppsurf_' strrep(filemot,'.mat','.fig')])
end

if isfield(F_map,'Pfe')
    surf(Id,Iq,Pfe), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('Iron Loss [W]')
    if not(kturns == 1)
        title(['Rewind factor = ' num2str(kturns)])
    end
    saveas(gcf,[NewDir 'Pfe' strrep(filemot,'.mat','.fig')])
end
