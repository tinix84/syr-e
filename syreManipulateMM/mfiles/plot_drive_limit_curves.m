% plot_drive_limit_curves

% torque
figure(10), hold on
plot(Plim.n,Plim.T,'k','LineWidth',[1.5])
plot(Plim.w_BC * rad2rpm,Plim.T_BC,'r','LineWidth',[1.5])
grid on, axis([0 Nmax 0 1.25*Tmax])
yy = ylabel('Nm'); xx = xlabel('rpm');
title('Torque')
% set(gca,'FontSize',[10],'FontWeight','bold')
h=gcf();
if isoctave()   %AS
    fig_name=strcat(pathname1, 'T');
    hgsave(h,[fig_name]);
else
    saveas(h,[pathname1 'T'])
end


% power
figure(11), hold on
plot(Plim.n,Plim.P,'k','LineWidth',[1.5])
plot(Plim.w_BC * rad2rpm,Plim.P_BC,'r','LineWidth',[1.5])
grid on, axis([0 Nmax 0 1.25*Pmax])
yy = ylabel('W'); xx = xlabel('rpm');
title('Shaft power')
h=gcf();
if isoctave()   %AS
    fig_name=strcat(pathname1, 'P');
    hgsave(h,[fig_name]);
else
    saveas(h,[pathname1 'P'])
end
%saveas(gcf,[pathname1 'P'])

% voltage
figure(12), hold on
plot(Plim.n,Plim.V*sqrt(3),'k','LineWidth',[1.5])
plot(Plim.w_BC * rad2rpm,Plim.V_BC*sqrt(3),'r','LineWidth',[1.5])
grid on, axis([0 Nmax 0 1.25*Vmax*sqrt(3)])
yy = ylabel('V'); xx = xlabel('rpm');
title('Line voltage - pk')
h=gcf();
if isoctave()   %AS
    fig_name=strcat(pathname1, 'V');
    hgsave(h,[fig_name]);
else
    saveas(h,[pathname1 'V'])
end
%saveas(gcf,[pathname1 'V'])

% current
figure(13), hold on, 
plot(Plim.n,Plim.I,'k','LineWidth',[1.5])
plot(Plim.w_BC * rad2rpm,Plim.I_BC,'r','LineWidth',[1.5]), grid on, %hold off
grid on, axis([0 Nmax 0 1.25*Imax])
yy = ylabel('A'); xx = xlabel('rpm');
title('Phase current - pk')
h=gcf();
if isoctave()   %AS
    fig_name=strcat(pathname1, 'I');
    hgsave(h,[fig_name]);
else
    saveas(h,[pathname1 'I'])
end
%saveas(gcf,[pathname1 'I'])

figure(14), subplot(2,1,1), hold on
plot(Plim.n,Plim.PF,'LineWidth',[1.5])
grid on, axis([0 Nmax 0.5 1])
yy = ylabel('PF'); xx = xlabel('rpm');
titolo = [motor_name ', Imax = ' num2str(Imax) 'A'];
titolo(titolo == '_') = '-';
title(titolo);
figure(14), subplot(2,1,2), hold on
plot(Plim.n,Plim.P,'k','LineWidth',[1.5]), hold on
plot(Plim.w_BC * rad2rpm,Plim.P_BC,'r','LineWidth',[1.5]), %hold off
grid on, axis([0 Nmax 0 1.25*Pmax])
yy = ylabel('Shaft Power [W]'); xx = xlabel('rpm');
legend(num2str(Imax_vect))
% title('Shaft power')
%adapt_figure_fonts('Times New Roman',12,10)h=gcf();
if isoctave()   %AS
    fig_name=strcat(pathname1, 'PF');
    hgsave(h,[fig_name]);
else
    saveas(h,[pathname1 'PF'])
end
%saveas(gcf,[pathname1 'PF'])

% print -dpsc2 AOA_CurrProf