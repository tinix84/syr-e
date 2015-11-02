% ...

clear all, close all

[FILENAME, PATHNAME, FILTERINDEX] = uigetfile('ImpProvaMTW*', 'CARICA DATI');
load([PATHNAME 'ImpProvaMTW']);
load([PATHNAME 'PuntiOtt']);
load([PATHNAME 'fdfq_idiq_n256']);

ultimo = min(DatiOpt.ultimo);

Pref = DatiOpt.Potenza(1,:);
nref = DatiOpt.velmec;
Ttab = DatiOpt.Tmap;


% coordinata potenza
iq_set = DatiOpt.IqMin(:,1:ultimo);
fd_set = interp2(Id,Iq,Fd,DatiOpt.IdMin(:,1:ultimo),DatiOpt.IqMin(:,1:ultimo));

Tref = [];
for jj = 1:length(nref)
    Tref(jj,:) = DatiOpt.Potenza(jj,1:ultimo) / (nref(jj)*pi/30);
    temp_iq_tab = interp1(Tref(jj,1:ultimo),iq_set(jj,1:ultimo),Ttab);
    temp_fd_tab = interp1(Tref(jj,1:ultimo),fd_set(jj,1:ultimo),Ttab);
    
    [Y,I] = max(temp_iq_tab);
    
    %% prolungamento curve - derivata ultimi due punti
    temp = diff(temp_iq_tab(1:I));
%     temp = temp(isfinite(temp));
    der_iq = mean(temp(end-2:end));
    temp_iq_tab(I+1:end) = temp_iq_tab(I) + (1:length(Tref)-I) * der_iq;
    temp_fd_tab(I+1:end) = temp_fd_tab(I);

%     temp_iq_tab(I+1:end) = temp_iq_tab(I);
%     temp_fd_tab(I+1:end) = temp_fd_tab(I);
    iq_tab(jj,:) = temp_iq_tab;
    fd_tab(jj,:) = temp_fd_tab;
end

iq_tab(iq_tab > 1050 * sqrt(2)) = 1050 * sqrt(2);

figure(1), hold off
p = plot(Ttab,fd_tab' * 5143 / 2,'-x'), grid on,
xlabel('Nm'), ylabel('fd\_set');
legend(p,num2str(nref(1)),num2str(nref(2)),num2str(nref(3)),num2str(nref(4)),num2str(nref(5)),num2str(nref(6)),num2str(nref(7)),num2str(nref(8)),num2str(nref(9)));

figure(2), hold off
p = plot(Ttab,iq_tab' /sqrt(2)*10,'-x'), grid on,
xlabel('Nm'), ylabel('iq\_set');
legend(p,num2str(nref(1)),num2str(nref(2)),num2str(nref(3)),num2str(nref(4)),num2str(nref(5)),num2str(nref(6)),num2str(nref(7)),num2str(nref(8)),num2str(nref(9)));

figure(3), hold off
p = plot(Ttab,fd_tab','-x'), grid on,
xlabel('Nm'), ylabel('fd\_set');
legend(p,num2str(nref(1)),num2str(nref(2)),num2str(nref(3)),num2str(nref(4)),num2str(nref(5)),num2str(nref(6)),num2str(nref(7)),num2str(nref(8)),num2str(nref(9)));

figure(4), hold off
p = plot(Ttab,iq_tab','-x'), grid on,
xlabel('Nm'), ylabel('iq\_set');
legend(p,num2str(nref(1)),num2str(nref(2)),num2str(nref(3)),num2str(nref(4)),num2str(nref(5)),num2str(nref(6)),num2str(nref(7)),num2str(nref(8)),num2str(nref(9)));

% stampa LUT
Tstep = diff(Ttab(1:2));
nstep = diff(nref(1:2));

% stampa
fid = fopen([PATHNAME 'fd_tab_iq_tab.txt'],'w');
fprintf(fid,['// IPM560 - ' date '\n']);
fprintf(fid,' -- \n');
fprintf(fid,'// TABELLE DEFLUSSAGGIO \n');
fprintf(fid,'// Tref min tabella    =  0 Nm\n');
fprintf(fid,'// Tref max tabella    =  %4.1f Nm\n',Ttab(end));
fprintf(fid,'// passo Tref          =  %4.1f Nm\n',Tstep);
fprintf(fid,'// n min tabella       =  %4.1f rpm\n',nref(1));
fprintf(fid,'// n max tabella       =  %4.1f rpm\n',nref(end));
fprintf(fid,'// passo n             =  %4.1f rpm\n',nstep);
fprintf(fid,'// Scala flussi: Vs x 5143\n');
fprintf(fid,'// Scala correnti: Arms x 10\n');
fprintf(fid,' -- \n');

fprintf(fid,'int ld_ref_tab[%2.0f][%2.0f]    =  {\n',length(nref),length(Ttab));
for jj = 1:size(fd_tab,2)
    fprintf(fid,'{');
    for kk = 1:(size(fd_tab,1)-1)
        fprintf(fid,'%6.0f,',fd_tab(kk,jj) * 5143);
    end
    fprintf(fid,'%6.0f}\n',fd_tab(end,jj) * 5143);
end
fprintf(fid,'}; \n');

fprintf(fid,'int iq_ref_tab[%2.0f][%2.0f]    =  {\n',length(nref),length(Ttab));
for jj = 1:size(iq_tab,2)
    fprintf(fid,'{');
    for kk = 1:(size(iq_tab,1)-1)
        fprintf(fid,'%6.0f,',iq_tab(kk,jj)/sqrt(2)*10);
    end
    fprintf(fid,'%6.0f}\n',iq_tab(end,jj)/sqrt(2)*10);
end
fprintf(fid,'}; \n');

fclose(fid);