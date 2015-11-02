
% plot_all_surfaces

% input:
% Id Iq Fd Fq
% optionally: loss maps

% output
% - surf of all maps in the input file (Id Iq Fd Fq and maybe loss maps)
% - surf of T, amp(F), amp(I), IPF

%% Fd Fq Id Iq
figure(1)
subplot(2,2,1)
surf(Id,Iq,Fd); grid on, hold on
xlabel('i_d [A]'),ylabel('i_q [A]'),zlabel('\lambda_d [Vs]')
% zlim([-0.1 0.15])
subplot(2,2,3)
surf(Id,Iq,Fq); grid on, hold on
xlabel('i_d [A]'),ylabel('i_q [A]'),zlabel('\lambda_q [Vs]'),
% zlim([0 0.25])
subplot(2,2,2)
surf(Id,Iq,Id); grid on, hold on
xlabel('i_d [A]'),ylabel('i_q [A]'),zlabel('i_d [A]')
% zlim([-0.1 0.15])
subplot(2,2,4)
surf(Id,Iq,Iq); grid on, hold on
xlabel('i_d [A]'),ylabel('i_q [A]'),zlabel('i_q [A]')
% zlim([0 0.25])
adapt_figure_fonts('Times New Roman',14,10)

%% F,I,T,IPF
figure(2)
subplot(2,2,1)
surf(Id,Iq,FI); grid on, hold on
xlabel('i_d [A]'),ylabel('i_q [A]'),zlabel('\lambda [Vs]')
% zlim([-0.1 0.15])
subplot(2,2,3)
surf(Id,Iq,TI); grid on, hold on
xlabel('i_d [A]'),ylabel('i_q [A]'),zlabel('T [Nm]'),
% zlim([0 0.25])
subplot(2,2,2)
surf(Id,Iq,II); grid on, hold on
xlabel('i_d [A]'),ylabel('i_q [A]'),zlabel('|i| [A]')
% zlim([-0.1 0.15])
subplot(2,2,4)
surf(Id,Iq,abs(IPF)); grid on, hold on
xlabel('i_d [A]'),ylabel('i_q [A]'),zlabel('I.P.F.')
% zlim([0 0.25])
adapt_figure_fonts('Times New Roman',14,10)

figure,
surf(fd,fq,ID_dato_fd_fq), grid on, xlabel('\lambda_d [Vs]'), ylabel('\lambda_q [Vs]')
zlabel('i_d [A]');% ylim([-0.1 0.25])
adapt_figure_fonts('Times New Roman',14,14)
axis([-0.1 0.1 0 0.2 -100 0])
saveas(gcf,[pathname 'Id surf']);

figure,
surf(fd,fq,IQ_dato_fd_fq), grid on, xlabel('\lambda_d [Vs]'), ylabel('\lambda_q [Vs]')
zlabel('i_q [A]');% ylim([-0.1 0.25])
adapt_figure_fonts('Times New Roman',14,14)
axis([-0.1 0.1 0 0.2 0 100])
saveas(gcf,[pathname 'Iq surf']);



