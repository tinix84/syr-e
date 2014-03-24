%% VALUTAZIONE DEI LIMITI SUPERIORI E INFERIORI PER LE VARIABILI DA OTTIMIZZARE

%  - valido per nlay qualunque
%  - il baricentro dei bounds è una macchina con:
%       * cave rotore equispaziate
%       * aria delle barriere, ponticelli e magneti distribuiti secondo una sinusoide in base alla funzione peso temp_sin_distribution
%  - il bound è nel range 0.5 - 2 rispetto al baricentro

% 23 04 10 - bounds_x

global tipo_cost

p = geo.p;                                      % Paia poli
nlay = geo.nlay;                                % Numero di barriere
ang_pont0 = geo.ang_pont0;                      % ampiezza (gradi) di un angolo che sottende un arco lungo pont0
pont0 = geo.pont0;                              % Ponticelli al traferro
x0 = geo.xr/geo.r;

% ANGOLI 'dalpha'

% Le (nlay+1) variabili chiamate 'dalpha' rappresentano gli angoli che definiscono la posizione delle barriere. Sono espressi in forma incrementale.  

temp_dalpha =  1/(nlay) * ones(1,nlay);
bounds_dalpha = temp_dalpha' * [1 2];
bounds_dalpha_1 = [16 (360/(4*p))/2*1.2];                    % bounds_dalpha_1 = [8 2*(360/(4*p))/(nlay+1)]; ????

% SPESSORI hc_pu ESPRESSI IN PER UNIT
temp_hc =  ones(1,nlay);
bounds_hc = temp_hc' * [0.2 1]; %14 aprile 0.2 invece di 0.1

% Induzione residua magneti 
    
bounds_Br = [0.2 0.8];

% Raggio al traferro (aggiunto 23 04 10)
bounds_x = x0 * [0.8  1.15];

% ANGOLO DI FASE DELLA CORRENTE ('gamma') 
    
bounds_gamma = [25 65];
% Numero di barriere di flusso ('nlayers') 
    
bounds_nlay = [0.51 4.49];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%temp_delta = 90/nlay/p; % passo regolare di dentatura (distribuzione uniforme su 90/p gradi, suddivisa in nlay)

%  1. dalpha
% distribuisco i dalpha aria come un seno e i ferri con quello che resta
% temp_dalpha_air = ones(1,nlay) * temp_delta;
% temp_dalpha_steel = ones(1,nlay) * temp_delta;
% temp_dalpha = zeros(1,length(temp_dalpha_steel) + length(temp_dalpha_air));
% temp_dalpha(1:2:end) = temp_dalpha_steel;
% temp_dalpha(2:2:end) = temp_dalpha_air;
% temp_dalpha = temp_delta * ones(1,2*geo.nlay) / geo.p;      % angolo
% meccanico
% temp_dalpha =  temp_delta * ones(1,2*geo.nlay);

% ang_min = geo.ang_pont0 * 2;
% bounds_dalpha(bounds_dalpha < ang_min) = ang_min;
% eccezione layer 1
% bounds_dalpha(1,:) = [10 45/2];
% somma dalpha
% bounds_dalpha_sum = [45-temp_delta*0.75 45-temp_delta*0.25];
% dalpha1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% PONTICELLI RADIALI -> IL LORO SPESSORE SI CALCOLA IN MODO DA RISPETTARE I VINCOLI DI RESISTENZA MECCANICA 

% last_pont = 1.0;
% bounds_pont = last_pont * ones(1,nlay-1)' * [0 1];
% bounds_pont = [bounds_pont;[last_pont last_pont]];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. pont
% temp_pont = temp_sin_distribution/temp_sin_distribution(1) * geo.pont0;
% bounds_pont = temp_pont' * [0.5 2];
% bounds_pont(bounds_pont > 2 * geo.pont0) = 2 * geo.pont0;

% bounds_pont = [0 1.0
%                0 1.0
%                1.0 1.0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% PERCENTUALE DI RIEMPIMENTO MAGNETE -> TUTTE LE BARRIERE VENGONO INTERAMENTE RIEMPITE DI MAGNETE

% bounds_magpu = ones(1,nlay)' * [0 1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. magpu
% temp_magpu = temp_sin_distribution/temp_sin_distribution(1) * 0.2;
% bounds_magpu = temp_magpu' * [0.5 2];
% bounds_magpu(bounds_magpu<0.1) = 0.1;
% bounds_magpu(bounds_magpu>1) = 1;
% bounds_magpu =[
%     0.1    1.0
%     0.1    1.0000
%     1.0    1.0000]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




