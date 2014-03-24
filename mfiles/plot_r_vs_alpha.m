% plot_r_vs_alpha

xr = 1; x0 = sqrt(2) * xr;
alpha = 0:1:45;

beta = 180/pi * calc_apertura_cerchio(pi/180*alpha,xr,x0);
r = (x0 - xr * cos(alpha*pi/180))./(cos(beta*pi/180)); 

figure(1)
plot(alpha,r), grid on,% hold on
% plot(alpha,alpha,'r'), hold off
legend('r/xr','alpha')
xlabel('\alpha - deg'), ylabel('r - p.u.'),