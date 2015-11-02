%% CAMBIO DI NOTAZIONE:NOTI I 2nlay ANGOLI 'dalpha',CALCOLO GLI nlay ANGOLI 'alpha' E LE nlay ALTEZZE 'hc',MEMORIZZANDOLE INOLTRE NELLA STRUTTURA 'geo'

% Input: dalpha
% Output: alpha, hc

% La nuova formulazione degli input geometrici prevede infatti di utilizzare gli angoli dalpha (sono gli angoli di posizione delle barriere, in
% forma incrementale) invece di far riferimento a alpha e hc separatamente.

function geo = calc_alpha_hc_delta_x0_2(geo)


xr = geo.xr;            % Raggio del rotore al traferro
p = geo.p;              % Paia poli
nlay = geo.nlay;        % N° layers
R = geo.r;              % Raggio ext
g = geo.g;              % Traferro
lt = geo.lt;            % Lunghezza denti
pont0 = geo.pont0;      % Ponticelli al traferro

%% MEMORIZZO DALLA STRUTTURA 'geo' GLI ANGOLI 'dalpha' E GLI SPESSORI 'hc_pu'
dalpha = geo.dalpha;
hc_pu = geo.hc_pu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NO:VEDI NUOVA NOTAZIONE IN PU+OPERAZIONE DI CONTROLLO SVOLTA DIRETTAMENTE IN OTTIMO
% se angolo totale eccessivo, riscalo tutto rispetto a un angolo
% sorteggiato (messo in bounds come ultimo dalpha)
% if sum(dalpha) > dalpha_tot
%     temp = dalpha_tot / sum(dalpha);
%     dalpha = temp * dalpha;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CALCOLO 'alpha' E LI SALVO IN GEO.
alpha = integra_fx(1:length(dalpha),dalpha);
geo.alpha = alpha;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure, plot(alpha_sum), hold on, plot(dalpha,'r'), hold off, keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CALCOLO LA POSIZIONE DEL CENTRO FITTIZIO x0, UTILE PER IL DISEGNO.
%% VALUTO IL RAGGIO DELL'ALBERO Ar.
x0 = xr/cos(pi/2/p);                            % centro cerchi (p generico)
geo.x0 = x0;
% 2013/07/25 GM tentativo tenendo fisso il valore dell'albero, la sua
% modifica porta di fatto a cabiare forma alla geo fluida.
Ar = x0 - xr * tan(pi/2/p);                     % Raggio dell'albero. (Il max valore consentito a questo parametro è calcolato in modo che il bordo
% geo.Ar = Ar;  %(decommentare per default)     % superiore dell'ultima barriera possa avere inizio, al limite, in corrispondenza della linea di
% terminazione del polo). Una volta calcolato, sovrascrivo Ar nella struttra geo.
htot = xr - Ar;                                 % spazio lineare disponibile

ly = R - xr - g - lt;                           % giogo statore
lyr = 1.0 * ly;                                 % ferro rotore garantito (tot)
la = xr - Ar -lyr;                              % massima aria disponibile

%% CALCOLO E MEMORIZZO LE ALTEZZE hc

% VALUTAZIONE DEGLI SPESSORI hc

hfe_min = 2 * pont0;                                                                % Spessore minimo guida di flusso
hc_half_min = 0.5; %0.5*pont0; limite minimo di hc 1 millimetro 14 aprile                                                                % Spessore minimo barriera di flusso

beta = 180/pi * calc_apertura_cerchio(pi/180*alpha,xr,x0);
r = (x0 - xr * cos(alpha*pi/180))./(cos(beta*pi/180));

hc = [];

if (nlay==1)

    %% max hc according to alpha min (26 Jan 2011)
    hc_half_max1 = (alpha*pi/180/(1+alpha*pi/180)*(xr-pont0)); 
    % (needs division by 2 .. don't know why but it works)
    hc_half_max1 = hc_half_max1 * 0.5;
    
    %% max hc according to alpha max (27 Jan 2011)
    temp_alpha_hfemin = hfe_min/xr; % rad
    temp_alpha_hc_2 = pi/(2*p) - alpha*pi/180 - temp_alpha_hfemin;
    hc_half_max2 = (temp_alpha_hc_2/(1+temp_alpha_hc_2)*(xr-pont0)); 
    hc_half_max = min(hc_half_max1,hc_half_max2);
    hc_pu(1) = 1;
    hc(1) = hc_pu(1) * hc_half_max * 2;
    if hc(1)<2*hc_half_min
        hc(1)=2*hc_half_min;
    end
    
    %% added on 25 Jan 2011
    if hc(1)>2*hc_half_max
        hc(1)=2*hc_half_max;
    end
    
else
    for jj = nlay:-1:1                                                                  % Parto dall'ultimo (do la precedenza al layer più interno)
        
        if (jj == nlay)&(jj>1) %attenzione tentativo di modifica per funzionare con 1 layer
            
            hc_half_max = (min((x0-r(end)-Ar),(r(end)-r(end-1)-hc_half_min)) - hfe_min);
            hc(jj) = hc_pu(jj) * hc_half_max * 2;
            
            if hc(jj)<2*hc_half_min
                hc(jj)=2*hc_half_min;
            end
            
        else
            
            if jj == 1
                
                hc_half_max = min((r(jj+1)-r(jj)-0.5*hc(jj+1)-hfe_min),(xr-pont0-x0+r(jj)));
                hc(jj) = hc_pu(jj) * hc_half_max * 2;
                
                if hc(jj)<2*hc_half_min && hc_half_min<=(xr-pont0-x0+r(jj))
                    hc(jj)=2*hc_half_min;
                end
            else
                
                hc_half_max = min((r(jj+1)-r(jj)-0.5*hc(jj+1)-hfe_min),(r(jj)-r(jj-1)-hc_half_min-hfe_min));
                hc(jj) = hc_pu(jj) * hc_half_max * 2;
                
                if hc(jj)<2*hc_half_min
                    hc(jj)=2*hc_half_min;
                end
            end
            
        end
        
    end
end

% SALVO hc NELLA STRUTTURA 'geo'.

geo.hc = hc;



%% CALCOLO 'nlay_effettivo', 'delta' E 'nr' E LI SALVO IN GEO.

geo.delta = [alpha(1) diff(alpha) (90/p)-alpha(end)];

geo.nr = ceil(90 ./ geo.delta ) * 2;