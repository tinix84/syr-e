% 25 05 07

% input:    kvMax_FdFq.mat
% output:   kvMax_idiq.mat

% addpath Mappe_IPMScooter\Mat

% clear

% load kvMax_FdFq
% load IdIq_FdFq

load([pathname 'kvMax_FdFq']);
load([pathname 'IdIq_FdFq']);

id_KvMax = interp2(fd,fq,ID_dato_fd_fq,fd_KvMax,fq_KvMax);
iq_KvMax = interp2(fd,fq,IQ_dato_fd_fq,fd_KvMax,fq_KvMax);

% save('mat\kvMax_idiq','id_KvMax','iq_KvMax','F_KvMax','T_KvMax')
save([pathname 'kvMax_idiq'],'id_KvMax','iq_KvMax','F_KvMax','T_KvMax');

figure;
plot(id_KvMax,iq_KvMax,'kx','LineWidth',[1.5])
grid on, xlabel('id [Arms]'),ylabel('iq [Arms]')
title(['MTPV (KvMax) - Imax = ' num2str(Imax/sqrt(2),4) 'Arms'],'FontSize',[12],'FontWeight','bold');
axis equal
