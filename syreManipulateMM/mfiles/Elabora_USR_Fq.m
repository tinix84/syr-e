% Elabora_USR_Fq
% Input:    Fq_0_step_7.m
%           Fq_448_step_7.m
%           dati importati da file .xw3 -> salvati come .m senza le intestazioni
% Output:   dati per CurveEstreme.mat (ex CurveFd.mat)            

% 1 - caratt magn Fq_0 (0,Iq)
K_V_Wave3 = 2760 / 256 / (380 * sqrt(2) / 1024);   %V pk stella to Wave3 voltage (duty cycle)   

% 31 07 07 -    
%               passo 7 --> Iqmax = 448Arms, Id = 0.
%               dati importati da file .xw3 
%               -> salvato come .m senza le intestazioni

run([PATHNAME curve_file3]);  % carica i dati

% Id = 0;

Exp.Id0 = Id * sqrt(2) * ones(1,64);
Exp.Iq0 = (64 * step : -step : step) * sqrt(2);

% Exp.Vd0 = [-1395 -1381 -1368 -1353 -1336 -1320 -1306 -1291 -1274 -1259 -1244 ...
%            -1226 -1209 -1195]/ K_V_Wave3;

%   complemento a 2, 16bit          

Exp.Vd0 = -(2^16 - v(1:64,2)) / K_V_Wave3;
Exp.Vd0 = Exp.Vd0';
Exp.Vd0 = fliplr(Exp.Vd0); 


Exp.Iq0 = Exp.Iq0(1:length(Exp.Vd0));

Exp.Vq0 = zeros(size(Exp.Vd0)) / K_V_Wave3;

Exp.Fq0 = - Exp.Vd0 / w;

%               passo 7 --> Iqmax = 448Arms, Id = 448.
%               dati importati da file .xw3 
%               -> salvato come .m senza le intestazioni

run([PATHNAME curve_file4]);  % carica i dati

Id = step * 64;

Exp.Id1 = Id * sqrt(2) * ones(1,64);
Exp.Iq1 = (64 * step : -step : step) * sqrt(2);

Exp.Vd1 = -(2^16 - v(1:64,2)) / K_V_Wave3;
Exp.Vd1 = Exp.Vd1';
Exp.Vd1 = fliplr(Exp.Vd1); 

Exp.Iq1 = Exp.Iq1(1:length(Exp.Vd1));
Exp.Vq1 = zeros(size(Exp.Vd1)) / K_V_Wave3;

Exp.Fq1 = - Exp.Vd1 / w;
% 
% save Fd_

figure
hold on
pl0 = plot(Exp.Iq0,Exp.Fq0,'b-x');
pl1 = plot(Exp.Iq1,Exp.Fq1,'r-x');
legend([pl0 pl1],'I_d = 0','I_d = 448Arms')
grid on
xlim([0 Id * sqrt(2)])
ylim([0 2.6])
title('Modello magnetico - Flusso q')
xlabel('I_q - Apk')
ylabel('Vs')


% salva le due curve per costruire il modello fg

Iq      = [Exp.Iq0 0]';
Iq0     = flipud(Iq);

FluxQo  = [Exp.Fq0 0]';
FluxQo  = flipud(FluxQo);

FluxQmax= [Exp.Fq1 0]';
FluxQmax= flipud(FluxQmax);
