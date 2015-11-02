% 25 05 07

% input:    ktMax_FdFq.mat
% output:   ktMax_idiq.mat

% addpath Mappe_IPMScooter

% load ktMax_FdFq
% load IdIq_FdFq

load([pathname 'ktMax_FdFq']);
load([pathname 'IdIq_FdFq']);

id_KtMax = interp2(fd,fq,ID_dato_fd_fq,fd_KtMax,fq_KtMax);
iq_KtMax = interp2(fd,fq,IQ_dato_fd_fq,fd_KtMax,fq_KtMax);

% save('mat\ktMax_idiq','id_KtMax','iq_KtMax','T_KtMax')
save([pathname 'ktMax_idiq'],'id_KtMax','iq_KtMax','T_KtMax')

figure;
[c,h]=contour(Id,Iq,T,max(max(T))*(0.1:0.1:0.9)); clabel(c,h); hold on
[c,h]=contour(Id,Iq,I,max(max(I))*(0.1:0.1:0.9)); clabel(c,h);
plot(id_KtMax,iq_KtMax,'kx','LineWidth',[1.5])
grid on, xlabel('id [Arms]'),ylabel('iq [Arms]')
title(['MTPA (KtMax) - Imax = ' num2str(Imax/sqrt(2),4) 'Arms'],'FontSize',[12],'FontWeight','bold');
axis equal

