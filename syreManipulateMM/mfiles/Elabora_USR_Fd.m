% Elabora_USR_Fd
% Input:    Fd_0_step_7.m
%           Fd_448_step_7.m
%           dati importati da file .xw3 -> salvati come .m senza le intestazioni
% Output:   dati per CurveEstreme.mat (ex CurveFd.mat)            

% 1 - caratt magn Fd_0 (Id,0)
K_V_Wave3 = 2760 / 256 / (380 * sqrt(2) / 1024);   %V pk stella to Wave3 voltage (duty cycle)   

% 31 07 07 -    Inom 500Arms
%               passo 7 --> Idmax = 448Arms, Iq = 0
%               dati importati da file .xw3 
%               -> salvato come .m senza le intestazioni

run([PATHNAME curve_file1]);  % carica i dati

w = N * pi/30 * 3;

Iq = 0;

Exp.Id0 = (64 * step : -step : step) * sqrt(2);
Exp.Iq0 = Iq * sqrt(2) * ones(1,64);

% Exp.Vq0 =  [3927 3915 3903 3891 3877 3864 3850 3838 3822 3809 3796 3779 3765 3751 ...
%             3736 3718 3702 3687 3671 3652 3635 3619 3601 3581 3563 3545 3526 3504 ...
%             3484 3464 3440 3420 3398] / K_V_Wave3
    
Exp.Vq0 = v(65:end,2) / K_V_Wave3;
Exp.Vq0 = Exp.Vq0';
Exp.Vq0 = fliplr(Exp.Vq0); 

Exp.Id0 = Exp.Id0(1:length(Exp.Vq0));
Exp.Vd0 = zeros(size(Exp.Vq0)) / K_V_Wave3;

Exp.Fd0 = Exp.Vq0 / w;

% 31 07 07 -    passo 7 --> Idmax = 448Arms, Iq = 448Arms
%               aumentate Imax e Inom di SR e IM
%               dati importati da file .xw3 
%               -> salvato come .m senza le intestazioni

run([PATHNAME curve_file2]);  % carica i dati

Iq = step * 64;           % rms

Exp.Id1 = (64 * step : -step : step) * sqrt(2);
Exp.Iq1 = Iq * sqrt(2) * ones(1,64);

% Exp.Vq1 =  [3808 3796 3783 3769 3752 3738 3724 3710 3691 ...
%             3677] / K_V_Wave3

Exp.Vq1 = v(65:end,2) / K_V_Wave3;
Exp.Vq1 = Exp.Vq1';
Exp.Vq1 = fliplr(Exp.Vq1); 

Exp.Id1 = Exp.Id1(1:length(Exp.Vq1));
Exp.Vd1 = zeros(size(Exp.Vq1)) / K_V_Wave3;

Exp.Fd1 = Exp.Vq1 / w;

figure
hold on
% pl0 = plot(Exp.Id0,Exp.Fd0,'b-x');
% pl1 = plot(Exp.Id1,Exp.Fd1,'r-x');
pl2 = plot(Exp.Id0,Exp.Fd0,'ro');
pl3 = plot(Exp.Id1,Exp.Fd1,'rx');

% 
legend([pl2 pl3],'I_q = 0',['I_q = ' num2str(Iq) 'Arms'])


grid on
% ylim([0 2.5]);
title('Modello magnetico - Flusso d');
xlabel('I_d - Apk');
ylabel('Vs');
xlim([0 Iq * sqrt(2)])
ylim([0 2.6])

% figure
% hold on
% plot(Exp.Id0/sqrt(2),Exp.Vq0/sqrt(2),'ro');
% plot(Exp.Id3/sqrt(2),-Exp.Vq1/sqrt(2),'bo');
% grid on
% xlabel('I_d - Arms');
% ylabel('V rms');


% salva le due curve per costruire il modello fg

Id      = [Exp.Id0 0]';
Id0     = flipud(Id);

FluxDo  = [Exp.Fd0 0]';
FluxDo  = flipud(FluxDo);

FluxDmax= [Exp.Fd1 0]';
FluxDmax  = flipud(FluxDmax);
