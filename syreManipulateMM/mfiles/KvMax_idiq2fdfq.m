% 08 11 09

% input:    kvMax_idiq.mat
% output:   kvMax_FdFq.mat

load([pathname 'kvMax_idiq']);

fd_KvMax = interp2(id,iq,Fd,id_KvMax,iq_KvMax);
fq_KvMax = interp2(id,iq,Fq,id_KvMax,iq_KvMax);

% save('mat\ktMax_idiq','id_KtMax','iq_KtMax','T_KtMax')
save([pathname 'kvMax_FdFq'],'fd_KvMax','fq_KvMax','T_KvMax')

% figure;
% [c,h]=contour(fd,fq,T,max(max(T))*(0.1:0.1:0.9)); clabel(c,h); hold on
% [c,h]=contour(fd,fq,I,max(max(I))*(0.1:0.1:0.9)); clabel(c,h);
% plot(fd_KvMax,fq_KvMax,'kx','LineWidth',[1.5])
% grid on, xlabel('fd [Vs]'),ylabel('fq [Vs]')
% title(['MTPV (KvMax) - Imax = ' num2str(Imax/sqrt(2),4) 'Arms'],'FontSize',[12],'FontWeight','bold');
% axis equal

