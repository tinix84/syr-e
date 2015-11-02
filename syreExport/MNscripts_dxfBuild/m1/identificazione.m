% identificazione.m
% - simula n punti di lavoro con corrente di ampiezza io per valutare il
% modello magnetico della macchina

function identificazione(geo,io,per);

% due posizioni per ciascuna corrente
nsim = 13;

passo_i = io / 5;
Id = -io:passo_i:io;
Iq = 0:passo_i:io;

risultati = [];

% curva lungo q, id = 0
gamma = 0; 
for jj = 1:length(Iq)
    modulo = Iq(jj);
    SOL = simula_60deg(geo,nsim,modulo,gamma);
    % risultati
    ris_sim = mean(SOL);
    T = abs(ris_sim(end));
    ris_sim(1) = gamma;
    risultati = [risultati; ris_sim];
    save identificazione risultati
end

size_q = jj;

% curva lungo d
gamma = -90;
for jj = 1:length(Id)
    modulo = Id(jj);
    SOL = simula_60deg(geo,nsim,modulo,gamma);
    % risultati
    ris_sim = mean(SOL);
    T = abs(ris_sim(end));
    ris_sim(1) = gamma;
    risultati = [risultati; ris_sim];
    save identificazione risultati
end

size_d = jj;

id = risultati(size_q+1:end,2);
iq = risultati(1:size_q,3);
fd = risultati(size_q+1:end,4);
fq = risultati(1:size_q,5);
% T = abs(risultati(:,end));

% f = sqrt(fd.^2 + fq.^2);
% delta = atan(fq./fd) * 180/pi;
% fi = delta - gamma;
% PF = cos(fi * pi/180);

% grafici
figure(1),
plot(id,fd,'-x'), grid on, hold on
plot(iq,fq,'-xr'), hold off
xlabel('id, iq - A'), ylabel('Flusso - Vs');
legend('\lambda_d','\lambda_q')
title('IDENTIFICAZIONE');

save identificazione risultati io id fd iq fq
