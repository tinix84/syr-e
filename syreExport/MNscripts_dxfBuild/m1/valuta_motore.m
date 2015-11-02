%% valuta_motore.m - 18 gen 10
% tratto da valuta_Tmaxnmax_simple.m

% lancia le simulazioni sulla base di geo e gamma_in

function [out] = valuta_motore(geo,per,io,gamma_in)

global tipo_valutazione options_global

if strcmp(tipo_valutazione,'gambj')
    gamma = gamma_in;
    nsim = geo.nsim_gambj; delta_sim = geo.delta_sim_gambj;
elseif strcmp(tipo_valutazione,'evalx')
    nsim = geo.nsim_evalx; delta_sim = geo.delta_sim_evalx;gamma = gamma_in;
elseif strcmp(tipo_valutazione,'singt')
    nsim = geo.nsim_singt; delta_sim = geo.delta_sim_singt;gamma = gamma_in;
elseif strcmp(tipo_valutazione,'singp') %%%%%%%%%123%%%%%%%%%%%%
    nsim = geo.nsim_singp; delta_sim = geo.delta_sim_singp; gamma = [0:20:40 45:5:65 70:10:90];
end

SOL = [];
T = zeros(1,length(gamma)); fd = T; fq = T; ripple = T;
Pfe_S_tot = T; Pfe_R_tot = T;

for kk = 1:length(gamma)
    SOL(:,:,kk) = simula_xdeg(geo,nsim,delta_sim,io,gamma(kk)); % risultati
    
    ris_sim = SOL(1:end-1,:,kk);
    
    T(kk) = abs(mean(ris_sim(:,6)));
    ripple(kk) = std(ris_sim(:,6));
    fd(kk) = mean(ris_sim(:,4));
    fq(kk) = mean(ris_sim(:,5));
   
  %% 22/01/2013 MG Commento il calcolo delle perdite nel ferro, allo stato attuale la presente geometria non è in grado di eseguire il calcolo correttamente
  %%
%     calc_Pfe1;
%     %keyboard
%     Pfe_S_tot(kk) = 2 * p * sum(Pfe_S);
%     Pfe_R_tot(kk) = 2 * p * sum(Pfe_R);
      Pfe_S_tot(kk)=0;
      Pfe_R_tot(kk)=0;
end

out.SOL = SOL(1:end-1,:,:);  %%ris_sim;
out.Tn = T;
out.ripple_pu = ripple./T;
out.fd = fd;
out.fq = fq;
out.Pfe_S = Pfe_S_tot;
out.Pfe_R = Pfe_R_tot;

% calcolo la potenza nominale (20 giugno 2010)
out.Pn=T/(abs(fd+1j*fq))*30*per.Vdc/(pi*sqrt(3)*geo.p);
% PER LA SIMULAZIONE A 3 OBIETTIVI DECOMMENTARE
    % Singola simulazione a 90° per valutare se fd=ld id
    SOL90 = SOL(end,:,1);
%    gamma = 90;
%    SOL90 = simula_xdeg(geo,2,delta_sim,io,gamma); % risultati
    ris_sim90 = (SOL90);
    T90 = abs(mean(ris_sim90(:,6)));
    fd90 = (mean(ris_sim90(:,4)));
    fq90 = (mean(ris_sim90(:,5)));

    out.SOL90 = SOL90;
    out.T90 = T90;
    out.fd90 = fd90;
    out.fq90 = fq90;
%     %modifica introdotta il 6 settembre 2010 da fc per calcolare il flusso dei
%     %magneti: si fa una ulteriore simulazione femm senza corrente
%     SOL_0A = simula_xdeg(geo,2,delta_sim,0,gamma); % risultati
%     ris_sim_0A = (SOL_0A);
%     fd_0A = (mean(ris_sim_0A(:,4)));
%     fq_0A = (mean(ris_sim_0A(:,5)));
% 
%     out.SOL_0A = SOL_0A;
%     out.fd_0A = fd_0A;
%     out.fq_0A = fq_0A;
%     out.fmagneti=sqrt(fd_0A.^2+fq_0A.^2);

%     %modifica introdotta il 13 marzo 2011 da fc per ????
%     gamma = 80;
%     SOL_80 = simula_xdeg(geo,2,delta_sim,io,gamma); % risultati
%     ris_sim_80 = (SOL_80);
%     T80 = abs(mean(ris_sim_80(:,end)));
%     fd_80 = (mean(ris_sim_80(:,4)));
%     fq_80 = (mean(ris_sim_80(:,5)));
%     out.SOL_80 = SOL_80;
%     out.T80 = T80;
%     out.fd_80 = fd_80;
%     out.fq_80 = fq_80;
%     out.fmagneti=sqrt(fd_80.^2+fq_80.^2);

% % % PER LA SIMULAZIONE A 3 OBIETTIVI COMMENTARE
% out.SOL90 = 0;
% out.T90 = 0;
% out.fd90 = 0;
% out.fq90 = 0;

% non usati
out.nn = 0;
out.nmax_pu = 0;

save sim_mot_temp out

