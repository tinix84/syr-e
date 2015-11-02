% 04 07 07

% Max torque profile

% input:    all the model (run C_AOA-...)
% output:   idq reference profile
%           Fdq reference profile

Fd_KtMax = interp2(id,iq,Fd,id_KtMax,iq_KtMax);
Fd_AB = interp2(id,iq,Fd,id_AB,iq_AB);
Fd_BC = interp2(id,iq,Fd,id_BC,iq_BC);

Fd_ref = [Fd_AB Fd_BC];

Fq_KtMax = interp2(id,iq,Fq,id_KtMax,iq_KtMax);
Fq_AB = interp2(id,iq,Fq,id_AB,iq_AB);
Fq_BC = interp2(id,iq,Fq,id_BC,iq_BC);

Fq_ref = [Fq_AB Fq_BC];

id_ref = [id_AB id_BC];
iq_ref = [iq_AB iq_BC];

w_ref = [w_AB w_BC] * rad2rpm;

[p_w_ref,s] = polyfit(1:length(w_ref),w_ref,5);
% w_ref_model=fit((1:length(w_ref))',w_ref','cubicspline')
w_ref = polyval(p_w_ref,1:length(w_ref));

Nn = ceil(w_AB(1) * rad2rpm / 10)*10;
% interpolo per LUT 16 punti (+1)
w_LUT = linspace(Nn,Nmax,17);
Fd_max = interp1(w_ref(1:end-6),Fd_ref(1:end-6),w_LUT);
iq_max = interp1(w_ref(1:end-6),iq_ref(1:end-6),w_LUT);
iq_max = iq_max / sqrt(2);


figure
hold on
plot(w_ref,Fd_ref,'LineWidth',[1.5]);
% plot(w_ref,Fq_ref,'r','LineWidth',[1.5]);
plot(w_LUT,Fd_max,'o','LineWidth',[1.5]);
grid on
xlim([0 Nmax]);
title('\lambda * - max power')
xlabel('rpm')
ylabel('Vs')

figure  % pu
hold on
plot(w_ref,Fd_ref * 5146.8,'LineWidth',[1.5]);
% plot(w_ref,Fq_ref,'r','LineWidth',[1.5]);
plot(w_LUT,Fd_max * 5146.8,'o','LineWidth',[1.5]);
grid on
xlim([0 Nmax]);
title('\lambda * - max power')
xlabel('rpm')
ylabel('Vs')

figure
hold on
plot(w_ref,iq_ref/sqrt(2),'LineWidth',[1.5]);
plot(w_LUT,iq_max,'o','LineWidth',[1.5]);
grid on
xlim([0 Nmax]);
title('i_q * - max power')
xlabel('rpm')
ylabel('Arms')

figure
hold on
plot(w_ref,sqrt(Fd_ref.^2 + Fq_ref.^2));
grid on
xlim([0 Nmax]);
title('Flux Amplitude - max power')
xlabel('rpm')
ylabel('Vs')


% save mat\LUT w_LUT iq_max Fd_max
save([PATHNAME 'LUT'],'w_LUT','iq_max','Fd_max')
% converto in valori numerici di macchina

% iq_max:   in decimi di Arms (x10)
% Fd_max:   KscalaF = 5146.8
% scala coppia

n_iq_max = iq_max * 10;
n_Fd_max = Fd_max * 5146.8;

StampaLUT;
% save mat\LUT_numerico w_LUT n_iq_max n_Fd_max
save([PATHNAME 'LUT_numerico'],'w_LUT','n_iq_max','n_Fd_max')
open([PATHNAME 'LUT.txt']);

