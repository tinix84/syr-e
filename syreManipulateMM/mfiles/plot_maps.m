
% plot_maps

% input:
% Id Iq Fd Fq
% optionally: loss maps

% output
% - mesh of all maps in the input file (Id Iq Fd Fq and maybe loss maps)
% - mesh of T, amp(F), amp(I), IPF

function plot_maps(pathname,filename)

if (nargin<1)
    pathname = cd;
    filename = 'fdfq_idiq_n256.mat'
end

load([pathname '\' filename])

% run([pathname 'ReadParameters']);

% Fd Fq curves
figure, hold all
[value, index] = min(abs(Iq(:,1)));
plot(Id(index,:),Fd(index,:)),
plot(Id(1,:),Fd(1,:)),
plot(Id(end,:),Fd(end,:)), hold on
[value, index] = min(abs(Id(1,:)));
plot(Iq(:,index),Fq(:,index)), hold all
plot(Iq(:,1),Fq(:,1)),
plot(Iq(:,end),Fq(:,end)),
xlabel('[A]'), ylabel('[Vs]')
% adapt_figure_fonts('Times New Roman',14,12)
saveas(gcf,[pathname 'fdfq curves'])

% Fd [Vs]
figure
mesh(Id,Iq,Fd), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('\lambda_d [Vs]')
hold on
plot3(Id(1,:),Iq(1,:),Fd(1,:),'k','LineWidth',2),
plot3(Id(end,:),Iq(end,:),Fd(end,:),'k--','LineWidth',2),
[value, index] = min(abs(Iq(:,1)));
plot3(Id(index,:),Iq(index,:),Fd(index,:),'k','LineWidth',2),
% adapt_figure_fonts('Times New Roman',14,12)
saveas(gcf,[pathname 'Fd mesh'])

% Fq [Vs]
figure
mesh(Id,Iq,Fq), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('\lambda_q [Vs]')
hold on
plot3(Id(:,1),Iq(:,1),Fq(:,1),'k','LineWidth',2),
plot3(Id(:,end),Iq(:,end),Fq(:,end),'k--','LineWidth',2),
[value, index] = min(abs(Id(1,:)));
plot3(Id(:,index),Iq(:,index),Fq(:,index),'k','LineWidth',2),
% if not(Kr == 1)
%     title(['Riavvolto ' num2str(Kr)])
% end
% adapt_figure_fonts('Times New Roman',14,12)
saveas(gcf,[pathname 'Fq mesh'])

% TORQUE
if exist('T','var')
    figure
    mesh(Id,Iq,abs(T)), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('Torque [Nm]')
%    adapt_figure_fonts('Times New Roman',14,12)
    saveas(gcf,[pathname 'Torque mesh'])
end
if exist('dT','var')
    figure
    mesh(Id,Iq,dT), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('Torque ripple [Nm]')
%     axis([0 30 0 30 0 0.5])
%    adapt_figure_fonts('Times New Roman',14,12)
    saveas(gcf,[pathname 'Torque Ripple mesh'])
end

% core loss
if exist('Pfes_h','var')
    %% Loss Map
    figure
    subplot(2,2,1)
    mesh(Id,Iq,Pfes_c); grid on, hold on
    xlabel('i_d [A]'),ylabel('i_q [A]'),zlabel('Pc-stat [W]')
    % zlim([-0.1 0.15])
    subplot(2,2,3)
    mesh(Id,Iq,Pfes_h); grid on, hold on
    xlabel('i_d [A]'),ylabel('i_q [A]'),zlabel('Ph-stat [W]'),
    % zlim([0 0.25])
    subplot(2,2,2)
    mesh(Id,Iq,Pfer_c); grid on, hold on
    xlabel('i_d [A]'),ylabel('i_q [A]'),zlabel('Pc-rot [W]')
    % zlim([-0.1 0.15])
    subplot(2,2,4)
    mesh(Id,Iq,Pfer_h); grid on, hold on
    xlabel('i_d [A]'),ylabel('i_q [A]'),zlabel('Ph-rot [W]')
    % zlim([0 0.25])
%    adapt_figure_fonts('Times New Roman',14,10)
    saveas(gcf,[pathname 'Loss mesh'])

end

% pm loss
if exist('Ppm','var')
    figure
    mesh(Id,Iq,Ppm), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('pm loss [W]')
%    adapt_figure_fonts('Times New Roman',14,12)
    saveas(gcf,[pathname 'PM loss mesh'])
end





