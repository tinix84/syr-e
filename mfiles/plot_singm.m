
% interp the flux linkage maps over a very dense grid
id = F_map.Id; iq = F_map.Iq;
Fd0 = F_map.Fd; Fq0 = F_map.Fq;
T = F_map.T;
dT = F_map.dT;

n = 256;
i_d=linspace(id(1),id(end),n);
i_q=linspace(id(1),iq(end),n);

[Id,Iq]=meshgrid(i_d,i_q);

Fd = interp2(id,iq,Fd0,Id,Iq,'cubic');
Fq = interp2(id,iq,Fq0,Id,Iq,'cubic');

T = interp2(id,iq,T,Id,Iq,'cubic');
T = T*klength;

dT = interp2(id,iq,dT,Id,Iq,'cubic');
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

save ([NewDir 'fdfq_idiq_n' num2str(n) '.mat'],'Fd','Fq','Id','Iq');
save ([NewDir 'fdfq_idiq_n' num2str(n) '.mat'],'T','-append');
save ([NewDir 'fdfq_idiq_n' num2str(n) '.mat'],'dT','-append');

% flux maps
figure
plot(Id(1,:),Fd([1 end],:)), grid on, hold on
plot(Iq(:,1),Fq(:,[1 end])), 
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

