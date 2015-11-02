% 12 01 08

% input:    kPFMax_FdFq.mat
% output:   kPFMax_idiq.mat

% addpath Mappe_IPMScooter\Mat

% clear

% load kvMax_FdFq
% load IdIq_FdFq

load([PATHNAME 'kPFMax_FdFq']);
load([PATHNAME 'IdIq_FdFq']);

id_PFMax = interp2(fd,fq,ID_dato_fd_fq,fd_PFMax,fq_PFMax);
iq_PFMax = interp2(fd,fq,IQ_dato_fd_fq,fd_PFMax,fq_PFMax);

% save('mat\kvMax_idiq','id_KvMax','iq_KvMax','F_KvMax','T_KvMax')
save([PATHNAME 'kPFMax_idiq'],'id_PFMax','iq_PFMax','PF_Max','T_PFMax');


figure
hold on
plot(id_KtMax,iq_KtMax,'b');
% plot(id_PFMax,iq_PFMax,'r');
grid on

