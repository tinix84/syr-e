

function DatiOpt = PuntiDiOttimo(d,setup,Id,Iq,Fd,Fq,Pfes_c,Pfer_c,Pfes_h,Pfer_h,Pmag,coeff_Pfe)

Vbus=d.Vbus;
Rs=d.Rs;
p=d.p;
velmin=d.velmin;
velbase=d.velbase;
velmax=d.velmax;
FreqEletDim=d.FreqEletDim;
PfeDim=d.PfeDim;
FlussoDim=d.FlussoDim;
Vbatt=d.Vbatt;
conf=d.conf;
velDim = d.velDim;

Kr = eval(setup{1});
Kl = eval(setup{2});
Kt = eval(setup{3});
nt = eval(setup{4});
ns = eval(setup{5});
motorType = (setup{6});
PMlossAfterSegmentation = eval(setup{7});

% PMlossAfterSegmentation=1;
% PMlossAfterSegmentation=0.7;
% PMlossAfterSegmentation=0.4;

%%%%%%%%%%%%%%%%%%
% change N turns, Kr = N/N0
Id=Id/Kr; Iq=Iq/Kr;
Fd=Fd*Kr; Fq=Fq*Kr;
d.Rs=d.Rs*Kr^2;

%%%%%%%%%%%%%%%%%%
% change stack length, Kl = stack/stack0
Fd=Fd*Kl;
Fq=Fq*Kl;
d.Rs = d.Rs*(Kl+Kt)/(1+Kt);

if exist('Pfes_c','var')
    Pfes_c = Pfes_c*Kl;
    Pfer_c = Pfer_c*Kl;
    Pfes_h = Pfes_h*Kl;
    Pfer_h = Pfer_h*Kl;
    Pmag = Pmag*Kl;
end

IDtot=Id(1,:);    % necessario per compatibilita'
IQtot=Iq(:,1)';   % necessario per compatibilita'

% torque [Nm]
T=3/2*p*(Fd.*Iq-Fq.*Id);

% torque range
Tmax=d.TMAX;

velmec = linspace(velmin,velmax,ns);
Tmap    = linspace(0,Tmax,nt);  % [Nm]
Tmap(1) = 0.05;                 % evitare il calcolo della coppia nulla

n_coppie=max(size(Tmap));
n_freq=max(size(velmec));

% elementi per interpolare ...
% n_flusso=10;
% n_Itau=10;

% Imposta il calcolo dell'ottimo in potenza sopra la velocita' base.
wbase = velbase*2*pi/60;          % [rad/s]
Pmap  = wbase*Tmap; %*Tbase/Tmap(end);    % [W]

% ps=4; % passo di stampa delle matrici finali
% pp=4; % passo di stampa delle matrici intermedie


% current vector
I=Id+1j*Iq;
% flusx linkage vector
F=Fd+1j*Fq;
% flux linkage amplitude
Fo=abs(F);

% Converter (not used)
% IniConvertitore;

FreqElet = p/60*velmec;                % [Hz]
wmecc=FreqElet*2*pi/p;                 % [rad/s]
rpm2rads=pi/30;
rads2Hz=1/(2*pi);
rpm2Hz=rpm2rads*rads2Hz;

% Torque contour curves
Curve=contourc(IDtot,IQtot',T,Tmap);


IdMin=[]; IqMin=[];
VoMin=[]; CosfiMin=[];
T_top_W=[];
ultimo=[];
P_min=[]; Potenza=[];
% PerdCONV=[];

PERDITE_JOULE=[];
PERDITE_FERRO_E_MAG=[];
PERDITE_FERRO=[];
if exist('Pfes_c','var')
    PERDITE_FERRO_SC=[];
    PERDITE_FERRO_SH=[];
    PERDITE_FERRO_SE=[];
    PERDITE_FERRO_RC=[];
    PERDITE_FERRO_RH=[];
    PERDITE_FERRO_RE=[];
end
PERDITE_MAG=[];
if exist('P_BarRot','var')
    PERDITE_BARRE_ROTORE=[];
end

% per ogni velocita'
for n = 1:n_freq,
    
    % 1. scale Pfe with frequency
    if exist('Pfes_c','var')
        % stator
        Pfe_cs = Pfes_c * (velmec(n)/velDim).^2;
        Pfe_hs = Pfes_h * (velmec(n)/velDim).^coeff_Pfe.alpha;
        % rotor
        Pfe_cr = Pfer_c * (velmec(n)/velDim).^2;
        Pfe_hr = Pfer_h * (velmec(n)/velDim).^coeff_Pfe.alpha;
        % stator + rotor
        Pfe = Pfe_cs + Pfe_hs + Pfe_cr + Pfe_hr;
    else
        % if loss map is not available: calc Pfe extrapolating a single
        % test value PfeDim at
        Pfe=PfeDim*(Fo/FlussoDim).^2*(FreqElet(n)/FreqEletDim)^1.8;
    end
    
    % 2. scale magnet loss (Perdite_mag) with frequency
    if exist('Pmag','var')
        if exist('coeff_Pmag','var')
            Perdite_mag = Pmag*interp1(coeff_Pmag.f_norm,coeff_Pmag.P_norm,(velmec(n)/velDim));
        else
            Perdite_mag=Pmag * (velmec(n)/velDim).^2;
        end
        Perdite_mag=PMlossAfterSegmentation*Perdite_mag;
    else
        Perdite_mag=zeros(size(Pfe));
    end
    
    % 3. back emf
    Vind=1j*2*pi*FreqElet(n)*F;
    
    % 4. current component representing Fe and PM loss (Ife_m)
    if (FreqElet(n)>0),
        Ife_m=2*(Pfe+Perdite_mag)./(3*abs(Vind));
        Ife=Ife_m.*exp(1j*angle(Vind));
    else
        Ife=0;
    end
    
    % aggiunge la corrente correlata alle perdite nel ferro alla corrente in ingresso
    Io=I+Ife;
    
    % 5. Motor voltage Voc = back emf + RI
    % Voc > Vbus is penalized by augmented loss, so to be excluded 
    Vof=Vind+Rs.*Io;    % phase voltage
    Voc=sqrt(3)*Vof;    % line voltage
    cosfi=cos(angle(Io)-angle(Vof));    % PF
    
    % calcola i moduli
    Voc_m=abs(Voc);
    Io_m= abs(Io);
    
    % total loss (3/2 goes with peak current values)
    Perdite=Pfe+Perdite_mag+3/2*Rs*Io_m.^2;
    if exist('P_BarRot', 'var')
        Perdite = Perdite+P_BarRot;
    end
      
    % (lim) equals 1 if Voc_m>Vbus and 0 otherwise
    lim=(sign(Voc_m-Vbus)+1)/2;
    lim(floor(Voc_m-Vbus)==0)=0;
    
    % augment loss if (lim)=(1)
    PerditeNew=Perdite+lim*1e10;
    
    % inizializzazione variabili ausiliarie
    Valori=0; Tot=0; idmin=[];     iqmin=[];
    vomin=[];     cosfimin=[];
    Pmin=[];     POT=[];
    PCONV=[];     fin=0;
    perdite_JOULE=[];
    perdite_FERRO_E_MAG=[];
    perdite_FERRO=[];
    if exist('Pfes_c','var')
        perdite_FERRO_SC=[]; perdite_FERRO_SH=[]; perdite_FERRO_SE=[];
        perdite_FERRO_RC=[]; perdite_FERRO_RH=[]; perdite_FERRO_RE=[];
    end
    perdite_MAG=[];
    if exist('P_BarRot','var')
        perdite_BARRE_ROTORE=[];
    end
    
    if (wmecc(n)>wbase)
        Curve=contourc(IDtot,IQtot',T*wmecc(n),Pmap);
    end
    
    % 6. torque (power) cycle:
    %    for each torque (power) setpoin, finds the id,iq condition that minimizes total loss  
    m = 0;
    while ((m<n_coppie)&&(Tot+Valori+1<length(Curve)))
        
        m = m+1;
        Tot=Tot+Valori+1;
        Valori=Curve(2,Tot);
        disp(Tot)
        disp(Valori)
        disp(' ')
        idIso=Curve(1,Tot+(1:Valori));
        iqIso=Curve(2,Tot+(1:Valori));
        PerditeIso=interp2(Id,Iq,PerditeNew,idIso,iqIso);
        
        % calcola la coppia di correnti che porta alla minima perdita
        [PminIso, indice]=min(PerditeIso);
        % assegna le correnti al percorso ottimo solo se il limite
        % di tensione e' rispettato
        
        if (PminIso<1e7)&&(fin<1),
            % (id,iq)-> Corrente esterna (solo nel caso frequenza==0, questa coincide con la corrente utile)
            id=idIso(indice);
            iq=iqIso(indice);
            
            %% debug
            ido = interp2(Id,Iq,real(Io),id,iq,'spline');
            iqo = interp2(Id,Iq,imag(Io),id,iq,'spline');
            if and(n>1,strcmp(motorType,'SR'))
                id = abs(ido);
            else
                id = ido;
            end
            iq = iqo;
            
            % variabili di ingresso costanti Inizializzate in IniConvertitore
            % Vbatt = tensione di batteria
            % conf  = configurazione del convertitore
            % componente = tipologia e caratteristiche del componente di potenza
            % condensatori = tipologia e caratteristiche dei condensatori
            
            VoP=interp2(Id,Iq,Voc_m,id,iq,'spline');    % Modulo della tensione concatenata
            IoP=interp2(Id,Iq,Io_m,id,iq,'spline');     % Modulo della corrente di linea
            cosfiP=interp2(Id,Iq,cosfi,id,iq,'spline'); % Fattore di potenza
            
            
            %             [Pcomm,Pcond,Pc_conv,Ib,IcrmsBATT,IcrmsBUS] = Perdconv(Vbatt,VoP,IoP,cosfiP,conf,componente,condensatori);
            %             PCONV=[PCONV Pcomm+Pcond+Pc_conv];
            
            if (wmecc(n)>=wbase)
                POT=[POT Pmap(m)];
            else
                POT=[POT wmecc(n)*Tmap(m)];
            end
            
            Pmin=[Pmin PminIso];
            idmin=[idmin id ];
            iqmin=[iqmin iq ];
            if (wmecc(n)>wbase)
                T_top=Pmap(m)/wmecc(n);
            else
                T_top=Tmap(m);
            end
            ult=m;
            
            if iqIso(indice)>max(max(Iq)),
                fin=1;
            end
            cosfimin=[cosfimin cosfiP];
            vomin=[vomin VoP];
            % 09-11-2010
            % Aggiungo mappo anche, per ciascun punto ottimo (velocità di
            % rotazione, coppia), i contributi di perdita: Joule, Ferro e
            % Mag
            provv_perdite = interp2(Id,Iq,(3*1/2*Rs*Io_m.^2),id,iq,'spline');
            perdite_JOULE = [perdite_JOULE provv_perdite];
            provv_perdite = interp2(Id,Iq,Pfe+Perdite_mag,id,iq,'spline');
            perdite_FERRO_E_MAG=[perdite_FERRO_E_MAG provv_perdite];
            provv_perdite = interp2(Id,Iq,Pfe,id,iq,'spline');
            perdite_FERRO=[perdite_FERRO provv_perdite];
            if exist('Pfes_c','var')
                provv_perdite = interp2(Id,Iq,Pfe_cs,id,iq,'spline');
                perdite_FERRO_SC=[perdite_FERRO_SC provv_perdite];
                provv_perdite = interp2(Id,Iq,Pfe_hs,id,iq,'spline');
                perdite_FERRO_SH=[perdite_FERRO_SH provv_perdite];
                %                 provv_perdite = interp2(Id,Iq,Pfe_es,id,iq,'spline');
                %                 perdite_FERRO_SE=[perdite_FERRO_SE provv_perdite];
                provv_perdite = interp2(Id,Iq,Pfe_cr,id,iq,'spline');
                perdite_FERRO_RC=[perdite_FERRO_RC provv_perdite];
                provv_perdite = interp2(Id,Iq,Pfe_hr,id,iq,'spline');
                perdite_FERRO_RH=[perdite_FERRO_RH provv_perdite];
                %                 provv_perdite = interp2(Id,Iq,Pfe_er,id,iq,'spline');
                %                 perdite_FERRO_RE=[perdite_FERRO_RE provv_perdite];
            end
            provv_perdite = interp2(Id,Iq,Perdite_mag,id,iq,'spline');
            perdite_MAG=[perdite_MAG provv_perdite];
            if exist('P_BarRot','var')
                provv_perdite = interp2(Id,Iq,P_BarRot,id,iq,'spline');
                perdite_BARRE_ROTORE=[perdite_BARRE_ROTORE provv_perdite];
            end
            
        else
            cosfimin=[cosfimin NaN];
            vomin=[vomin NaN];
            PCONV=[PCONV NaN];
            POT=[POT NaN];
            Pmin=[Pmin NaN];
            idmin=[idmin NaN];
            iqmin=[iqmin NaN];
            perdite_JOULE=[perdite_JOULE NaN];
            perdite_FERRO_E_MAG=[perdite_FERRO_E_MAG NaN];
            perdite_FERRO=[perdite_FERRO NaN];
            if exist('Pfes_c','var')
                perdite_FERRO_SC=[perdite_FERRO_SC NaN];
                perdite_FERRO_SH=[perdite_FERRO_SH NaN];
                %                 perdite_FERRO_SE=[perdite_FERRO_SE NaN];
                perdite_FERRO_RC=[perdite_FERRO_RC NaN];
                perdite_FERRO_RH=[perdite_FERRO_RH NaN];
                %                 perdite_FERRO_RE=[perdite_FERRO_RE NaN];
            end
            perdite_MAG=[perdite_MAG NaN];
            if exist('P_BarRot','var')
                perdite_BARRE_ROTORE = [perdite_BARRE_ROTORE NaN];
            end
        end
    end
    ultimo(n)=ult;   % contiene l'indice dell'ultimo elemento utile prima dei NaN
    T_top_W=[T_top_W; T_top]; % contiene la massima coppia realizzabile
    IdMin=[IdMin; idmin]; % percorso ottimo di Id
    IqMin=[IqMin; iqmin]; % percorso ottimo di Iq
    P_min=[P_min; Pmin];  % perdite minime
    CosfiMin= [CosfiMin; cosfimin];
    VoMin=[VoMin; vomin];
    Potenza=[Potenza; POT];
    %     PerdCONV=[PerdCONV; PCONV];
    PERDITE_JOULE=[PERDITE_JOULE; perdite_JOULE];
    PERDITE_FERRO_E_MAG=[PERDITE_FERRO_E_MAG; perdite_FERRO_E_MAG];
    PERDITE_FERRO=[PERDITE_FERRO; perdite_FERRO];
    if exist('Pfes_c','var')
        PERDITE_FERRO_SC=[PERDITE_FERRO_SC; perdite_FERRO_SC];
        PERDITE_FERRO_SH=[PERDITE_FERRO_SH; perdite_FERRO_SH];
        %         PERDITE_FERRO_SE=[PERDITE_FERRO_SE; perdite_FERRO_SE];
        PERDITE_FERRO_RC=[PERDITE_FERRO_RC; perdite_FERRO_RC];
        PERDITE_FERRO_RH=[PERDITE_FERRO_RH; perdite_FERRO_RH];
        %         PERDITE_FERRO_RE=[PERDITE_FERRO_RE; perdite_FERRO_RE];
    end
    PERDITE_MAG=[PERDITE_MAG; perdite_MAG];
    if exist('P_BarRot','var')
        PERDITE_BARRE_ROTORE=[PERDITE_BARRE_ROTORE; perdite_BARRE_ROTORE];
    end
    
end


DatiOpt.ultimo   = ultimo;
DatiOpt.T_top_W  = T_top_W;
DatiOpt.IdMin    = IdMin;
DatiOpt.IqMin    = IqMin;
DatiOpt.P_min    = P_min;
DatiOpt.CosfiMin = CosfiMin;
DatiOpt.VoMin    = VoMin;
DatiOpt.Potenza  = Potenza;
% DatiOpt.PerdCONV = PerdCONV;
DatiOpt.PERDITE_JOULE=PERDITE_JOULE;
DatiOpt.PERDITE_FERRO_E_MAG=PERDITE_FERRO_E_MAG;
DatiOpt.PERDITE_FERRO=PERDITE_FERRO;
if exist('Pfes_c','var')
    DatiOpt.PERDITE_FERRO_SC=PERDITE_FERRO_SC;
    DatiOpt.PERDITE_FERRO_SH=PERDITE_FERRO_SH;
    DatiOpt.PERDITE_FERRO_SE=PERDITE_FERRO_SE;
    DatiOpt.PERDITE_FERRO_RC=PERDITE_FERRO_RC;
    DatiOpt.PERDITE_FERRO_RH=PERDITE_FERRO_RH;
    DatiOpt.PERDITE_FERRO_RE=PERDITE_FERRO_RE;
end
DatiOpt.PERDITE_MAG=PERDITE_MAG;
if exist('P_BarRot','var')
    DatiOpt.PERDITE_BARRE_ROTORE=PERDITE_BARRE_ROTORE;
else
    DatiOpt.PERDITE_BARRE_ROTORE=zeros(size(PERDITE_JOULE));
end

DatiOpt.Tmap=Tmap;
DatiOpt.velmec=velmec;
% DatiOpt.componente=componente;
% DatiOpt.condensatori=condensatori;
DatiOpt.IDtot=IDtot;
DatiOpt.IQtot=IQtot;

DatiOpt.Kr=Kr;          % Riavvolgimento
DatiOpt.Kl=Kl;          % Allungamento

% elaborazioni dei risultati
if (1)
    % MOTOR
    DatiOpt.PotCONV = DatiOpt.Potenza + DatiOpt.P_min;
    DatiOpt.EffMOT  = 100 * DatiOpt.Potenza./DatiOpt.PotCONV;
    %   DatiOpt.PerdSYS = DatiOpt.P_min + DatiOpt.PerdCONV;
    %   DatiOpt.EffCONV = 100 * PotCONV./(PotCONV+ DatiOpt.PerdCONV);
    %   DatiOpt.EffSYS  = 0.01*EffCONV.*EffMOT;
else
    % GENERATOR
    DatiOpt.EffMOT  = 100 * (DatiOpt.Potenza - DatiOpt.P_min)./DatiOpt.Potenza;
    DatiOpt.PotCONV = DatiOpt.Potenza - DatiOpt.P_min;
    DatiOpt.PerdSYS = DatiOpt.P_min + DatiOpt.PerdCONV;
    DatiOpt.EffCONV = 100 * (PotCONV - DatiOpt.PerdCONV)./PotCONV;
    DatiOpt.EffSYS  = 0.01*EffCONV.*EffMOT;
end


