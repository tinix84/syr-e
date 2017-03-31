

function Mac=genMacDataFromGeo(geo,mat)

Mac.g=geo.g;
Mac.R=geo.R;
Mac.RtS=geo.r+geo.g;
Mac.p=geo.p;
Mac.Q=geo.ns*geo.p;
Mac.q=Mac.Q/(6*geo.p);
Mac.l=geo.l;
Mac.Qs=geo.Qs;
Mac.ps=geo.ps;
Mac.th0=geo.th0;
Mac.N_turn=geo.Ns;
Mac.N_cond=geo.Ns/geo.p/Mac.q/size(geo.avv,1); % numero di conduttori in cava per strato
Mac.N_parallel=1;
Mac.avv=geo.avv;
Br = mat.LayerMag.Br;
if (Br==0.0)
    Mac.n_mag_simulati=0;
else
    Mac.n_mag_simulati=size(geo.BLKLABELS.rotore.BarName,1);
end
if strcmp(geo.RotType,'SPM')
    Mac.n_mag_simulati=geo.ps;
end
Mac.fem=geo.fem;
Mac.caso=1;
