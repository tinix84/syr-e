
% identificazione.m
% - simula n punti di lavoro con corrente di ampiezza io per valutare il
% modello magnetico della macchina

function F_map = eval_FdFq_tables_in_FEMM(geo,id_vett,iq_vett,CurveEstreme)

nsim = round(geo.nsim_singt/2);
delta_sim = geo.delta_sim_singt;

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

risultati = [];

I = Id + 1i * Iq;
io = abs(I);
io_femm=io*geo.Nbob;
gamma = angle(I) * 180/pi;

Fd = zeros(size(Id)); Fq = Fd; T = Fd; dT = Fd;

for rr = 1:size(io,1)
    for cc = 1:size(io,2)
        SOL = simulate_xdeg(geo,nsim,delta_sim,io_femm(rr,cc),gamma(rr,cc));
        dT(rr,cc) = std(SOL(1:end-1,6));
        ris_sim = mean(SOL(1:end-1,:),1);
        T(rr,cc) = (ris_sim(6));
        Id(rr,cc)=(ris_sim(2))/geo.Nbob;
        Id(rr,cc)=(ris_sim(3))/geo.Nbob;
        Fd(rr,cc) = (ris_sim(4))*geo.Nbob;
        Fq(rr,cc) = (ris_sim(5))*geo.Nbob;
    end
end

F_map.Id = Id;
F_map.Iq = Iq;
F_map.Fd = Fd;
F_map.Fq = Fq;
F_map.T = T;
F_map.dT = dT;

% save sim_mot_temp F_map