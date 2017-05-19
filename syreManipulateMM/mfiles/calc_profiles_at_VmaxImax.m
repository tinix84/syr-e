function Plim = calc_profiles_at_VmaxImax(pathname,FILENAME,imax,rot_temperature)

debug = 0;


load([pathname FILENAME]);
if exist([pathname 'ReadParameters.m'])
    run([pathname 'ReadParameters']);
end

rad2rpm=30/pi/p;

% load MTPA and MTPV files
pathname = [pathname,'AOA\'];
load([pathname 'ktMax_idiq']); load([pathname 'kvMax_idiq']);

% Performance maps
id = Id(1,:); iq = Iq(:,1)';
TI = 3/2 * p * (Fd .* Iq - Fq .* Id);   % torque
FI = sqrt(Fd.^2 + Fq.^2);               % flux linkage amplitude
II = sqrt(Id.^2 + Iq.^2);               % current amplitude
IPF = cos(atan(Fq./Fd)-atan(-Id./Iq));  % internal PF

FImax = max(max(abs(FI)));
% no load flux
lm = interp2(Id,Iq,FI,0,0,'cubic');
if isnan(lm)
    lm = 0;
end

if debug
    figure (100)
    hold on
    plot(id_KvMax,iq_KvMax,'-x');
    plot(id_KtMax,iq_KtMax,'r-x'), hold off
    grid on
end

if strcmp(axes_type,'SR')
    % poly fit KtMax
    [p_KtMax_i,s] = polyfit(id_KtMax,iq_KtMax,7);
    [p_KtMax_T,s] = polyfit(id_KtMax,T_KtMax,7);
    id_KtMax_p = linspace(0,max(id_KtMax),length(id_KtMax));
    iq_KtMax_p = polyval(p_KtMax_i,id_KtMax_p);
    T_KtMax_p = polyval(p_KtMax_T,id_KtMax_p);
    % poly fit KvMax
    if ~isempty(id_KvMax)
        [p_KvMax_i,s] = polyfit(id_KvMax,iq_KvMax,3);
        id_KvMax = linspace(0,1.10*max(id_KvMax),length(id_KvMax))
        iq_KvMax = polyval(p_KvMax_i,id_KvMax);
    end
else
    % poly fit KtMax
    [p_KtMax_i,s] = polyfit(iq_KtMax,id_KtMax,7);
    [p_KtMax_T,s] = polyfit(iq_KtMax,T_KtMax,7);
    iq_KtMax_p = linspace(0,max(iq_KtMax),length(iq_KtMax));
    id_KtMax_p = polyval(p_KtMax_i,iq_KtMax_p);
    T_KtMax_p = polyval(p_KtMax_T,iq_KtMax_p);
    % poly fit KvMax
    [p_KvMax_i,s] = polyfit(iq_KvMax,id_KvMax,3);
    iq_KvMax = linspace(0,1.15*max(iq_KvMax),length(iq_KvMax));
    id_KvMax = polyval(p_KvMax_i,iq_KvMax);
end

if (debug)
    figure (100) ,hold on
    plot(id_KvMax,iq_KvMax,'k');
    plot(id_KtMax,iq_KtMax,'k'), hold off, axis equal
    grid on
    pause
end

% definitions:
% point A - MTPA @ I = imax
% point B - start of MTPV, still @ I = imax
% punto C - characteristic current

I_KtMax = abs(id_KtMax + j* iq_KtMax);
I_KvMax = abs(id_KvMax + j*iq_KvMax);

if not(isempty(id_KvMax))
    ich = min(I_KvMax);
else
    ich = inf;
end

% point A
if strcmp(axes_type,'SR')
    
    id_A = interp1(I_KtMax,id_KtMax,imax);
    iq_A = polyval(p_KtMax_i,id_A);
    %iq_A = interp1(I_KtMax,iq_KtMax,imax);
    
  
    if (ich>imax)
        % no MTPV or no MTPV and Imax crossing
        id_B = 0;
        iq_B = imax;
        id_C = id_B;
        iq_C = iq_B;
    else % SR style axes
        id_B = interp1(I_KvMax,id_KvMax,imax);
        iq_B = polyval(p_KvMax_i,id_B);
        %iq_B = interp1(I_KvMax,iq_KvMax,imax);
        id_C = 0;
        iq_C = polyval(p_KvMax_i,id_C);
    end
    
else
    % PM style axes
    iq_A = interp1(I_KtMax,iq_KtMax,imax);
    id_A = polyval(p_KtMax_i,iq_A);
    %id_A = interp1(I_KtMax,id_KtMax,imax);
    
    if (ich>imax)
        % no MTPV or no MTPV and Imax crossing
        id_B = -imax;
        iq_B = 0;
        id_C = id_B;
        iq_C = iq_B;
    else
        iq_B = interp1(I_KvMax,iq_KvMax,imax);
        id_B = polyval(p_KvMax_i,iq_B);
        %id_B = interp1(I_KvMax,id_KvMax,imax);
        iq_C = 0;
        id_C = polyval(p_KvMax_i,iq_C);
    end
end


% F_A = interp2(id,iq,FI,id_A,iq_A,'spline',0);  % rated flux
% T_A = interp2(id,iq,TI,id_A,iq_A,'spline',0);  % rated torque

F_A = interp2(id,iq,FI,id_A,iq_A);  % rated flux
T_A = interp2(id,iq,TI,id_A,iq_A);  % rated torque


w1 = Vmax / F_A;                    % base speed - elt rad/s
w1_rpm = w1 * rad2rpm;

F_B = interp2(id,iq,FI,id_B,iq_B);  % b = 0, imax flux
T_B = interp2(id,iq,TI,id_B,iq_B);  % rated torque

w2 = Vmax / F_B;                    % elt rad/s
w2_rpm = w2 * rad2rpm

% tratto A --> B
c = contourc(id,iq,II,[imax imax]);

id_Imax = (c(1,2:end));
iq_Imax = (c(2,2:end));

if strcmp(axes_type,'SR')
    % A --> B
    id_AB = id_Imax((iq_Imax >= iq_A) & (iq_Imax <= iq_B));
    iq_AB = iq_Imax((iq_Imax >= iq_A) & (iq_Imax <= iq_B)); %avoid interp errors
    % B --> C
    if (ich<=imax)
        id_BC = linspace(id_B,id_C,50);
        id_BC = id_BC(1:end-1);
        %iq_BC = polyval(p_KvMax_i,id_BC);
        iq_BC = interp1(id_KvMax,iq_KvMax,id_BC);
    else
        id_BC = [id_B id_C];
        iq_BC = [iq_B iq_C];
    end
else
    % SPM style axes
    id_AB = id_Imax((iq_Imax >= iq_B) & (iq_Imax <= iq_A));
    iq_AB = iq_Imax((iq_Imax >= iq_B) & (iq_Imax <= iq_A)); %avoid interp errors
    % tratto B --> C
    if (ich<=imax)
        iq_BC = linspace(iq_B,iq_C,50);
        iq_BC = iq_BC(1:end-1);
        %id_BC = polyval(p_KvMax_i,iq_BC);
        iq_BC = interp1(iq_KvMax,id_KvMax,iq_BC);
    else
        id_BC = [id_B id_C];
        iq_BC = [iq_B iq_C];
    end
end

% clear id_Imax iq_Imax

% PROFILO COPPIA - POT MAX
F_AB = interp2(id,iq,FI,id_AB,iq_AB);

if F_AB(end) > F_AB(1)
    F_AB = fliplr(F_AB);
    id_AB = fliplr(id_AB);
    iq_AB = fliplr(iq_AB);
end
w_AB = Vmax ./ F_AB;
T_AB = interp2(id,iq,TI,id_AB,iq_AB);
V_AB = w_AB .* F_AB;
I_AB = interp2(id,iq,II,id_AB,iq_AB);

F_BC = interp2(id,iq,FI,id_BC,iq_BC);
w_BC = Vmax ./ F_BC;
T_BC = interp2(id,iq,TI,id_BC,iq_BC);
V_BC = w_BC .* F_BC;
I_BC = interp2(id,iq,II,id_BC,iq_BC);

% low speed values (w < w1)
w_0A = linspace(0,w1,20);
T_0A = ones(size(w_0A)) * T_A;
V_0A = F_A * w_0A;
I_0A = ones(size(w_0A)) * imax;

% limiti iq
id_max = [id_A*ones(size(w_0A)) id_AB id_BC];
iq_max = [iq_A*ones(size(w_0A)) iq_AB iq_BC];

F_max = [F_A*ones(size(w_0A)) F_AB F_BC];
iq_min = zeros(size(iq_max));

% a pieno carico
fd_AB = interp2(id,iq,Fd,id_AB,iq_AB);
fd_BC = interp2(id,iq,Fd,id_BC,iq_BC);
fq_AB = interp2(id,iq,Fq,id_AB,iq_AB);
fq_BC = interp2(id,iq,Fq,id_BC,iq_BC);
%  a vuoto
fq_AB_0 = interp2(id,iq,Fq,id_AB,0);
fq_BC_0 = interp2(id,iq,Fq,id_BC,0);
% fq_AB_0 = interp2(id,iq,Fd,id_AB,0);
% fq_BC_0 = interp2(id,iq,Fd,id_BC,0);
% a vuoto

fd_max = [fd_AB(1)*ones(size(w_0A)) fd_AB fd_BC];
fq_max = [fq_AB(1)*ones(size(w_0A)) fq_AB fq_BC];

fq_0 = [fq_AB_0(1)*ones(size(w_0A)) fq_AB_0 fq_BC_0];

F0 = abs(fd_max + j*fq_0);
temp = abs(fd_max + j*fq_max);
iq_min(lm > F_max) = iq_max(lm > F_max);

% mechanical speed evaluation (IM only)
% synchronous speed
if exist('Wslip','var')
    wslip = interp2(Id,Iq,Wslip,id_max,iq_max);
    [a,last_number] = find(not(isnan(wslip)),1,'last');
    wslip(isnan(wslip)) = wslip(last_number);
    %     rot_temperature = 20;
    %     ref_temperature = 100;
    temp_coeff = (234.5 + rot_temperature)/(234.5 + Rr_temp);
    wslip = wslip * temp_coeff;
else
    wslip = id_max * 0;
end

w = [w_0A w_AB w_BC];
wr = w - wslip;

% potenza meccanica corretta (a parte Pfe)
P = [T_0A T_AB T_BC] .* wr/p;
% calcolo sbagliato, manca wslip .. lascio per compatibilità con versioni prec
P_AB = T_AB .* w_AB/p;
P_BC = T_BC .* w_BC/p;
P_0A = T_0A .* w_0A/p;

% Tmax = T_A;
% Pmax = max([P_AB P_BC])

PFlim = P ./ (3/2 * [V_0A V_AB V_BC] .* [I_0A I_AB I_BC]);


% save Plim points
Plim.wr = wr;
Plim.n = wr * rad2rpm;
Plim.P = P;
Plim.V = [V_0A V_AB V_BC];
Plim.I = [I_0A I_AB I_BC];
Plim.T = [T_0A T_AB T_BC];
Plim.F = F_max;
Plim.fd_max = fd_max;
Plim.fq_max = fq_max;
Plim.id_max = id_max;
Plim.iq_max = iq_max;

Plim.PF = PFlim;

Plim.iq_min = iq_min;

Plim.id_A = id_A;
Plim.iq_A = iq_A;
Plim.id_B = id_B;
Plim.iq_B = iq_B;
Plim.F_A = F_A;
Plim.F_B = F_B;
Plim.T_A = T_A;

Plim.w_BC = w_BC;
Plim.I_BC = I_BC;
Plim.V_BC = V_BC;
Plim.P_BC = P_BC;
Plim.T_BC = T_BC;

