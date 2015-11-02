% print_model
% run by BuildCurveEstreme
% fig 

figure
hold on
pl1 = plot(Id/sqrt(2),FluxDo,'ro');
pl2 = plot(Id/sqrt(2),FluxDmax,'rx');
pl3 = plot(Iq/sqrt(2),FluxQo,'o');
pl4 = plot(Iq/sqrt(2),FluxQmax,'x');
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