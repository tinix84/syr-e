
id = F_map.Id; iq = F_map.Iq;
Fd0 = F_map.Fd; Fq0 = F_map.Fq;
T = F_map.T;
dT = F_map.dT;

n = 256;
i_d=linspace(id(1),id(end),n);
i_q=linspace(id(1),iq(end),n);

[Id,Iq]=meshgrid(i_d,i_q);

fd = interp2(id,iq,Fd0,Id,Iq,'cubic');
fq = interp2(id,iq,Fq0,Id,Iq,'cubic');

T = interp2(id,iq,T,Id,Iq,'cubic');
T = T*Klenght;

dT = interp2(id,iq,dT,Id,Iq,'cubic');
dT = dT*Klenght;

%% rewind
Id=Id/kturns;
Iq=Iq/kturns;
Ld=Ld*kturns;
Lq=Lq*kturns;
%% lunghezza
Ld=Ld*Klenght;
Lq=Lq*Klenght;

% %% add end-connections term
% Ld = Ld + Lld * Id;
% Lq = Lq + Llq * Iq;

    save ([pathname 'ldlq_idiq_n' num2str(PuntiQuadro) '.mat'],'Ld','Lq','Id','Iq');
    if exist('T','var')
        save ([pathname 'ldlq_idiq_n' num2str(PuntiQuadro) '.mat'],'T','-append');
    end
    if exist('dT','var')
        save ([pathname 'ldlq_idiq_n' num2str(PuntiQuadro) '.mat'],'dT','-append');
    end


figure
surfc(Id,Iq,Ld), grid on, xlabel('id'), ylabel('iq'), zlabel('\lambda_d')
if not(Kr == 1)
    title(['Riavvolto ' num2str(Kr)])
end
saveas(gcf,[pathname1 'Fdsurf_' motor_name '.fig'])

figure
surfc(Id,Iq,Lq), grid on, xlabel('id'), ylabel('iq'), zlabel('\lambda_q')
if not(Kr == 1)
    title(['Riavvolto ' num2str(Kr)])
end
saveas(gcf,[pathname1 'Fqsurf_' motor_name '.fig'])


%% TORQUE MAP
    figure
    surf(Id,Iq,abs(T)), grid on, xlabel('id [A]'), ylabel('iq [A]'), zlabel('Torque [Nm]')
    if not(Kr == 1)
        title(['Riavvolto ' num2str(Kr)])
    end
    saveas(gcf,[pathname1 'Torquesurf_' motor_name '.fig'])

    figure
    surf(Id,Iq,dT), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('Torque ripple [Nm]')
    if not(Kr == 1)
        title(['Riavvolto ' num2str(Kr)])
    end
    axis([0 30 0 30 0 0.5])
    saveas(gcf,[pathname1 'Torqueripplesurf_' motor_name '.fig'])

