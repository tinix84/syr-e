% visualizza_risultati
clear all, close all
load risultati

torque = SOL(:,4:end);
figure(1)
plot(th,torque,'-x'), grid on, hold on
plot(th,mean(torque'),'k'), hold off
xlabel('posizione angolare'),ylabel('Nm')
title('torque')
gcf

T = repeat_n(mean(torque'),6);
h = spettro_pu(T,60,2);

figure(3)
plot(T,'k'), hold off
xlabel('posizione angolare'),ylabel('Nm')
title('torque')
gcf

fdq = SOL(:,7:8)
f123 = SOL(:,4:6)
th = SOL(:,1)
figure(4)
plot(th,SOL(:,4:8),'x-'), grid on