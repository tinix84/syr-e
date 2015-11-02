% identificazione.m
% - simula n punti di lavoro con corrente di ampiezza io per valutare il
% modello magnetico della macchina

function identificazione1(geo,id_vett,iq_vett,CurveEstreme)

nsim = 16;
delta_sim = 30;

%% PM-like axes
id_vett = -id_vett;

idm = min(id_vett); iqm = min(iq_vett);
idM = max(id_vett); iqM = max(iq_vett);
if (CurveEstreme)
    temp1 = id_vett;
    temp2 = temp1;
    temp3 = idm*ones(size(temp1));
    temp4 = idM * ones(size(temp1));
    Id = [temp1;temp2;temp3;temp4];
    temp1 = iqm * ones(size(temp1));
    temp2 = iqM * ones(size(temp1));
    temp3 = iq_vett;
    temp4 = temp3;
    Iq = [temp1;temp2;temp3;temp4];
else
    temp1 = id_vett;
    temp2 = iq_vett;
    [Id,Iq] = meshgrid(temp1,temp2);
end
%keyboard
risultati = [];

I = Id + 1i * Iq;
io = abs(I);
gamma = angle(I) * 180/pi - 90;

Fd = zeros(size(Id)); Fq = Fd; T = Fd; dT = Fd;
%keyboard
for rr = 1:size(io,1)
    for cc = 1:size(io,2)
        %keyboard
        SOL = simula_xdeg(geo,nsim,delta_sim,io(rr,cc),gamma(rr,cc));
        %keyboard
        dT(rr,cc) = std(SOL(1:end-1,6));
        ris_sim = mean(SOL(1:end-1,:),1);
        T(rr,cc) = (ris_sim(6));
        Fd(rr,cc) = (ris_sim(4));
        Fq(rr,cc) = (ris_sim(5));
        if size(ris_sim,2)>6
            numero_barriere = 0.5*(size(ris_sim,2)-6);
            for bb=1:numero_barriere
                BLavoroPM(rr,cc,bb) = abs(ris_sim(6+bb)+1i*ris_sim(6+numero_barriere+bb));
            end
        end
%         save sim_mot_temp Id Iq Fd Fq T dT
    end
end
%keyboard
F_map.Id = Id;
F_map.Iq = Iq;
F_map.Fd = Fd;
F_map.Fq = Fq;
F_map.T = T;
F_map.dT = dT;
F_map.BLavoroPM = BLavoroPM;
save sim_mot_temp F_map