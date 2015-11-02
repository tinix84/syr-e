% 08 11 09

% input:    ktMax_idiq.mat
% output:   ktMax_FdFq.mat

load([pathname 'ktMax_idiq']);
% load([pathname 'IdIq_FdFq']);

fd_KtMax = interp2(id,iq,Fd,id_KtMax,iq_KtMax);
fq_KtMax = interp2(id,iq,Fq,id_KtMax,iq_KtMax);

% save('mat\ktMax_idiq','id_KtMax','iq_KtMax','T_KtMax')
save([pathname 'ktMax_FdFq'],'fd_KtMax','fq_KtMax','T_KtMax')

% figure;
% [c,h]=contour(fd,fq,T,max(max(T))*(0.1:0.1:0.9)); clabel(c,h); hold on
% [c,h]=contour(fd,fq,I,max(max(I))*(0.1:0.1:0.9)); clabel(c,h);
% plot(fd_KtMax,fq_KtMax,'kx','LineWidth',[1.5])
% grid on, xlabel('fd [Vs]'),ylabel('fq [Vs]')
% title(['MTPA (KtMax) - Imax = ' num2str(Imax/sqrt(2),4) 'Arms'],'FontSize',[12],'FontWeight','bold');
% axis equal

