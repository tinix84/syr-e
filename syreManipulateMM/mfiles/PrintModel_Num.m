% print_model
% run by BuildCurveEstreme
% fig 

K_Vs_Fref = 5.1468e+003;    % Vs to Hitachi Flux units

figure
hold on
pl1 = plot(Id/sqrt(2),FluxDo * K_Vs_Fref,'ro');
pl2 = plot(Id/sqrt(2),FluxDmax * K_Vs_Fref,'rx');
pl3 = plot(Iq/sqrt(2),FluxQo * K_Vs_Fref,'o');
pl4 = plot(Iq/sqrt(2),FluxQmax * K_Vs_Fref,'x');
legend([pl1 pl2 pl3 pl4],'I_q = 0',['I_q = ' ...
    num2str(max(Iq)/sqrt(2)) 'Arms'],'I_d = 0' ...
    ,['I_d = ' num2str(max(Id)/sqrt(2)) 'Arms']);
grid on
% ylim([0 2.5]);
title('Modello magnetico');
xlabel('I_d, I_q - Arms');
ylabel('Vs');
xlim([0 max(Iq)/sqrt(2)])
ylim([0 2.6])