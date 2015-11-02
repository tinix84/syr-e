% function [StatAveHyst,StatAveEddy,RotAveHyst,RotAveEddy] = (geo,out,per)

% Calcolo perdite nel ferro
%
% Parametri dell'equazione dei Steinmetz
% P = K_h * f^alpha * B^beta + K_e * (sfB)^2

% Carica documento - CANCELLARE
% clear all;close all; clc;
% pathName = ([cd '\ProvaIronLoss\mot_65_T_eval_13_60d\']);
% docName = ('mot_65_T_eval_13_60d.fem');
% solName = ('mot_65_T_eval_13_60d.ans');
% structName = ('mot_65_T_eval_13_60d.mat');
% openfemm;
% opendocument([pathName docName]);
% opendocument([pathName solName]);
% load([pathName structName]);
% nPtsLoss = 5;

% M530-65A
% alpha = 1.45253; 
% beta = 1.90274;
% kh = 0.00689662;
% ke = 3.46139e-005;

% MassDensity = 7700;                 % Mass density [kg/m^3] 

% JFE 10JNEX900
alpha = geo.alpha;
beta = geo.beta;
kh = geo.kh;
ke = geo.ke;
 
MassDensity = geo.rhoFE;                 % Mass density [kg/m^3]

Rst = geo.R;                        % Stator External Radius [mm]
Rist = geo.r + geo.g;               % Stator Internal Radius [mm]
lt = geo.lt;                        % Thooth length [mm]
wt = geo.wt;                        % Thooth width [mm]
Rrot = geo.r;                       % Rotor radius [mm]
Lax = geo.l;                        % Axial Length [mm]
nSl = geo.q*geo.p*2*3;              % Number of Stator Slots
wbase = geo.nmax;                       % Rotor speed [rpm]
p = geo.p*2;                        % Number of poles
ps = geo.ps;                        % Number of simulated poles
Qs = geo.Qs;                        % Number of simulated slots
f_el = (wbase*((2*pi)/60)*geo.p)/(2*pi);        % electrical frequency [Hz]
deltaSim = geo.delta_sim_singt;     %


nSim = length(out.FluxDens(1,:))-3;         % Simulation number
nn = length(out.FluxDens(:,1));             % Mesh elements number
coordMeshElem = out.FluxDens(:,1);          % Mesh elements centroid coordinates as complex number
areaMeshElem = out.FluxDens(:,2);           % Mesh elements area using the same length unit defined for the problem [mm^2]
groupMeshElem = out.FluxDens(:,3);          % Mesh element group (1 Stator - 2 Rotor)

% Recover the b info
bMeshElem = zeros(nn,nSim);
clear ii; clear jj;
for ii = 1:nn
    for jj = 1:nSim
        bMeshElem(ii,jj) = out.FluxDens(ii,jj+3);
    end
end
bMeshElem = bMeshElem';                     

% Rebuilt the entire electrical period
clear ii;
clear jj;
Bx_actual = zeros(nSim,nn);
By_actual = zeros(nSim,nn);
for ii = 1:nSim
    for jj = 1:nn
        Bx_actual(ii,jj) = real(bMeshElem(ii,jj));
        By_actual(ii,jj) = imag(bMeshElem(ii,jj));
    end
end

Bx_actualM = zeros(nSim+1,nn);
BxInv = zeros(nn,(nSim*360/deltaSim)+1);
By_actualM = zeros(nSim+1,nn);
ByInv = zeros(nn,(nSim*360/deltaSim)+1);

clear ii

Bx = zeros(nSim*2,nn);
By = zeros(nSim*2,nn);
for ii = 1:nn
    if(groupMeshElem(ii)) == 1
        Bx(:,ii) = [Bx_actual(:,ii);-Bx_actual(:,ii)];
        By(:,ii) = [By_actual(:,ii);-By_actual(:,ii)];
    elseif(groupMeshElem(ii)) == 2
        Bx(:,ii) = [Bx_actual(:,ii);Bx_actual(:,ii)];
        By(:,ii) = [By_actual(:,ii);By_actual(:,ii)];
    end
end

Bx = [Bx;Bx(1,:)];
By = [By;By(1,:)];

bTotMeshElem = zeros(size(Bx));
clear ii; clear jj;
for ii = 1:length(Bx(:,1));
    for jj = 1:nn
        bTotMeshElem(ii,jj) = Bx(ii,jj)+j*By(ii,jj);
    end
end

% FFT computed as amplitude of each harmonic at the centroid
% of each element in the mesh.

clear ii; clear jj;
for jj = 1:nn
    [BxFFT_ho(:,jj),BxFFT_mag(:,jj),BxFFT_phase(:,jj)] = FFTAnalysis(real(bTotMeshElem(:,jj)),50,0);  % faccio l'analisi fino alla 50^ armonica
    [ByFFT_ho(:,jj),ByFFT_mag(:,jj),ByFFT_phase(:,jj)] = FFTAnalysis(imag(bTotMeshElem(:,jj)),50,0);  % faccio l'analisi fino alla 50^ armonica
end
BFFT = sqrt((BxFFT_mag.^2)+(ByFFT_mag.^2));

volMeshElem = (areaMeshElem*Lax)*1e-9;                           % Volume [m^3]
                              
weightMeshElem = volMeshElem*MassDensity;                        %  Weight [kg]

clear ii; clear jj;
for ii = 1:length(BxFFT_ho(:,1))
    for jj = 1:nn
        Physt(ii,jj) = kh*((BxFFT_ho(ii,jj)*f_el)^alpha)*(BFFT(ii,jj)^beta);        % [W/kg]
        Peddy(ii,jj) = ke*((BxFFT_ho(ii,jj)*f_el)^2*BFFT(ii,jj)^2);                     % [W/kg]
    end
end

%%
% keyboard
theta = (0:360/(length(bTotMeshElem(:,1))-1):360);

clear ii; clear jj;
for ii = 1:nn
        deltaBx(ii) = BxFFT_mag(2,ii)*2;
        deltaBxpp(ii) = deltaBx(ii)^(beta-alpha);
        polBx(:,ii) = polyfit(theta',real(bTotMeshElem(:,ii)),3);
        derBx(:,ii) = polyder(abs(polBx(:,ii)));
        fBx(:,ii) = kh*derBx(:,ii).^alpha*deltaBxpp(ii)*wbase*pi/30;
        intBx(:,ii) = polyint(fBx(:,ii)');
        PhystSullx(:,ii) = (1/(2*pi))*polyval(intBx(:,ii),2*pi)-polyval(intBx(:,ii),0);
        
        deltaBy(ii) = ByFFT_mag(2,ii)*2;
        deltaBypp(ii) = deltaBy(ii)^(beta-alpha);
        polBy(:,ii) = polyfit(theta',imag(bTotMeshElem(:,ii)),3);
        derBy(:,ii) = polyder(abs(polBy(:,ii)));
        fBy(:,ii) = kh*derBy(:,ii).^alpha*deltaBypp(ii)*wbase*pi/30;
        intBy(:,ii) = polyint(fBy(ii)');
        PhystSully(:,ii) = (1/(2*pi))*polyval(intBy(:,ii),2*pi)-polyval(intBy(:,ii),0);
end

PhystSull = PhystSullx+PhystSully; 


clear ii; clear jj;
for ii = 1:length(Physt(:,1))
    for jj = 1:nn
        Physt_W(ii,jj) = weightMeshElem(jj)*Physt(ii,jj);                       % [W]
        Peddy_W(ii,jj) = weightMeshElem(jj)*Peddy(ii,jj);                       % [W]
    end
end

Tot_Physt_W = sum(Physt_W);                  % [W]
Tot_Peddy_W = sum(Peddy_W);                  % [W]
Tot_Physt_W = Tot_Physt_W';                  % [W]
Tot_Peddy_W = Tot_Peddy_W';                  % [W]
 
clear ii;
 for ii = 1:nn
    if groupMeshElem(ii) == 1
        PhystStat_W(ii) = Tot_Physt_W(ii);                 % [W]
        PeddyStat_W(ii) = Tot_Peddy_W(ii);                 % [W]
    else
        PhystStat_W(ii) = 0;                 % [W]
        PeddyStat_W(ii) = 0;                 % [W]
    end
 end
 
 clear ii;
 for ii = 1:nn
    if groupMeshElem(ii) == 2
        PhystRot_W(ii) = Tot_Physt_W(ii);                  % [W]
        PeddyRot_W(ii) = Tot_Peddy_W(ii);                  % [W]
    else
        PhystRot_W(ii) = 0;                  % [W]
        PeddyRot_W(ii) = 0;                  % [W]
    end
 end
 
StatTotPhyst = sum(PhystStat_W)*p/ps;
StatTotPeddy = sum(PeddyStat_W)*p/ps;
RotTotPhyst = sum(PhystRot_W)*p/ps;
RotTotPeddy = sum(PeddyRot_W)*p/ps;


IronLoss = [StatTotPhyst StatTotPeddy RotTotPhyst RotTotPeddy];         % [Stat_Hyst_W Stat_Eddy_W Rot_Hyst_w Rot_Eddy_W]

%% DA CANCELLARE - SOLO PER VERIFICA PERDITE TOTALI

 clear ii; clear jj;
 for ii = 1:length(BxFFT_ho(:,1))
 for jj = 1:nn
    if groupMeshElem(jj) == 2
        VerPhystRot_W(ii,jj) = Physt_W(ii,jj);                  % [W]
        VerPeddyRot_W(ii,jj) = Peddy_W(ii,jj);                  % [W]
    else
        VerPhystRot_W(ii,jj) = 0;                  % [W]
        VerPeddyRot_W(ii,jj) = 0;                  % [W]
    end
 end
 end

 VerTotPhystRot = sum((VerPhystRot_W*(p/ps))');
 VerTotPeddyRot = sum((VerPeddyRot_W*(p/ps))');
 VerTotPRot = VerTotPhystRot+VerTotPeddyRot;
 
  clear ii; clear jj;
 for ii = 1:length(BxFFT_ho(:,1))
 for jj = 1:nn
    if groupMeshElem(jj) == 1
        VerPhystStat_W(ii,jj) = Physt_W(ii,jj);                  % [W]
        VerPeddyStat_W(ii,jj) = Peddy_W(ii,jj);                  % [W]
    else
        VerPhystStat_W(ii,jj) = 0;                  % [W]
        VerPeddyStat_W(ii,jj) = 0;                  % [W]
    end
 end
 end

 VerTotPhystStat = sum((VerPhystStat_W*(p/ps))');
 VerTotPeddyStat = sum((VerPeddyStat_W*(p/ps))');
 VerTotPRot = VerTotPhystStat+VerTotPeddyStat;
 
%%

