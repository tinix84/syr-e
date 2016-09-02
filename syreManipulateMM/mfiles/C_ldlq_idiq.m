
%% incremental Fd Fq
dFdd = diff(Fd,1,2)./diff(Id,1,2);
dFdq = diff(Fd,1,1)./diff(Iq,1,1);

dFqq = diff(Fq,1,1)./diff(Iq,1,1);
dFqd = diff(Fq,1,2)./diff(Id,1,2);

mask_data = (abs(Id(1,:))>0);

errAngCross = atand(2*dFdq(:,2:end)./(dFdd(2:end,:)-dFqq(:,2:end)));
ksnslss = (dFdd(2:end,:).*dFqq(:,2:end) - dFqq(:,2:end).^2 - 2*dFdq(:,2:end).^2)./(dFdd(2:end,:).*dFqq(:,2:end) - dFdq(:,2:end).^2);

figure
surfc(Fd); grid on, hold on
figure
surfc(Id(2:end,2:end),Iq(2:end,2:end),dFdd(2:end,:)); grid on, hold on, title('dFd/id')
figure
surfc(Id(2:end,2:end),Iq(2:end,2:end),dFdq(:,2:end)); grid on, hold on, title('dFd/iq')

nds = 6; % downsampling
figure(1001), 
subplot(2,2,1)
plot(Id(1,1:end-1),dFdd(1:nds:end,:)), grid on, xlabel('i_d [A]'), ylabel('l_d [H]'); %ylim([0 0.004])
subplot(2,2,2)
plot(Id(1,1:end),dFdq(1:nds:end,:)), grid on, xlabel('i_d [A]'),ylabel('l_{dq} [H]'); %ylim([-0.002 0.002])
subplot(2,2,3)
plot(Iq(1:end,1),dFqd(:,1:nds:end)), grid on, xlabel('i_q [A]'),ylabel('l_{qd} [H]'); %ylim([-0.002 0.002])
subplot(2,2,4)
plot(Iq(1:end-1,1),dFqq(:,1:nds:end)), grid on, xlabel('i_q [A]'),ylabel('l_q [H]'); %ylim([0 0.004])
adapt_figure_fonts('Times New Roman',12,10)
saveas(gcf,[pathname 'Fdq_inc']);

% cross saturation error
figure
errAngCross(errAngCross>10) = 10;
errAngCross(errAngCross<-10) = -10;
[C,h] = contour(Id(1,2:end),Iq(2:end,1),errAngCross); grid on, hold on, title('error [elt deg]')
clabel(C,h)
xlabel('id [A]'), ylabel('iq [A]'),
axis([0 10 0 10]), axis equal

%% apparent dFdq, approx #1 

LLd = (Fd-fm)./Id;
LLq = Fq ./ Iq;

figure
[c,h] = contour(Id,Iq,LLd,[0.001 0.002 0.003 0.004 0.0045 0.005 0.0055 0.006]), grid on, axis equal, clabel(c,h)
xlabel('i_d [A]'), ylabel('i_q [A]'), axis([-35 0 0 35])
adapt_figure_fonts('Times New Roman',12,10)

figure
[c,h] = contour(Id,Iq,LLq), grid on, axis equal, clabel(c,h)
xlabel('i_d [A]'), ylabel('i_q [A]'), axis([-35 0 0 35])
adapt_figure_fonts('Times New Roman',12,10)


% mask_data = (abs(Id(1,:))>0);

% figure
% surf(LLd), hold on, surf(LLq)

figure(1011),% title('d inductance - apparent','FontSize',[12],'FontWeight','boFd');
subplot(2,2,1)
plot(Id(1,:),LLd(1:nds:end,:)), grid on, xlabel('i_d [A]'),ylabel('L_d [H]');
% ylim([0 0.004])
subplot(2,2,2)
plot(Id(1,:),Id(1,:)*0),
grid on, xlabel('i_d [A]'),ylabel('L_{dq} [H]');
% ylim([-0.002 0.002])
subplot(2,2,3)
plot(Iq(:,1),Iq(:,1)*0),
grid on, xlabel('i_q [A]'),ylabel('L_{qd} [H]');
% ylim([-0.002 0.002])
subplot(2,2,4)
plot(Iq(:,1),LLq(:,1:nds:end)), grid on, xlabel('i_q [A]'),ylabel('L_q [H]');
% ylim([0 0.004])
adapt_figure_fonts('Times New Roman',12,10)
saveas(gcf,[pathname 'Fdq_app1']);

% apparent dFdq, approx #2
LLLd = (Fd(1,:)-fm)./Id(1,:);
LLLd = ones(length(LLLd),1)*LLLd;
LLLdq = LLd - LLLd;

LLLq = Fq(:,1)./Iq(:,1);
LLLq = LLLq*ones(1,length(LLLq));
LLLqd = LLq - LLLq;

figure(1021),% title('d inductance - apparent','FontSize',[12],'FontWeight','boFd');
subplot(2,2,1)
plot(Id(1,:),LLLd(1:nds:end,:)), grid on, xlabel('i_d [A]'),ylabel('L_d [H]'); ylim([0 0.004])
subplot(2,2,2)
plot(Id(1,:),LLLdq(1:nds:end,:)); grid on, xlabel('i_d [A]'),ylabel('L_{dq} [H]'); ylim([-0.002 0.002])
subplot(2,2,3)
plot(Iq(:,1),LLLqd(:,1:nds:end)), grid on, xlabel('i_q [A]'),ylabel('L_{qd} [H]'); ylim([-0.002 0.002])
subplot(2,2,4)
plot(Iq(:,1),LLLq(:,1:nds:end)), grid on, xlabel('i_q [A]'),ylabel('L_q [H]'); ylim([0 0.004])
adapt_figure_fonts('Times New Roman',12,10)
saveas(gcf,[pathname 'Fdq_app2']);
