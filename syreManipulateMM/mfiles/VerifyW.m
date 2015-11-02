% verify_W
% input:    CurveEstreme.mat
% output:   Co-energy plot

% verifica coenergia EXP
% occhio ai nomi    - FluxQmax  è Id = 0
%                   - FluxQo    è Id = max
Fd0_model = fit(Id,FluxDo,'splineinterp');
Fdmax_model = fit(Id,FluxDmax,'splineinterp');
Fq0_model = fit(Iq,FluxQo,'splineinterp');
Fqmax_model = fit(Iq,FluxQmax,'splineinterp');

int_Fd0 = integrate(Fd0_model, Id, 0);
int_Fdmax = integrate(Fdmax_model, Id, 0);
int_Fq0 = integrate(Fq0_model, Iq, 0);
int_Fqmax = integrate(Fqmax_model, Iq, 0);

figure
pl = plot(Id / sqrt(2),[int_Fd0-int_Fdmax int_Fq0-int_Fqmax])
grid on
hold on
legend(pl,'int\_\lambda_{d0}','int\_\lambda_{dmax}','int\_\lambda_{q0}','int\_\lambda_{qmax}')
% plot(Iq_Exp,[int_Fq0 int_Fqmax])
xlabel('Arms')
text(500,500,['\DeltaW''_d = ' num2str(int_Fd0(end) - int_Fdmax(end))])
text(500,400,['\DeltaW''_q = ' num2str(int_Fqmax(end) - int_Fq0(end))])
title('EXP - Integrale delle Curve Estreme - Id,Iq = 0 -> 448 Arms.','FontSize',[12],'FontWeight','Bold')
xlim([0 900])