% calc_io.m
% data la geometria e la dissipazione ammessa calcola l'ampiezza del vettore corrente
% input:
% - Kj : W /m2 dissipabili sullo statore
% - g : geometria della macchina (+ car. avvolgimento)
% output:
% - io (A)

% ref. Tutorial Course Notes: Design, Analysis and Control of Interior PM Synchronous Machines
% IEEE IAS Annual Meeting 2004 - Seattle
% cap. 6 - Vagati

function io = calc_io(geo,per)

% qui dentro si usa notazione normalizzata stile Guglielmi - Vagati
% per cui non badare troppo a significato dei simboli

% superficie esterna statore - m2
Asup = 2*pi*geo.r*geo.l*1e-6; 
% dissipazione - W/m2
kj = per.Loss / Asup; 
% temperatura
tempcu = per.tempcu;
% paia poli
p = geo.p;
% induzione normalizzata nel ferro
b = geo.b;
% coeff. scaling dente
kt = geo.kt;

% spire in serie per fase
n_spire = geo.Ns;
% raccorciamento di passo
kracc = geo.kracc;
% riempimento cava
kcu = geo.kcu;

% raggio statore NORMr
r = 1;
% lunghezza pacco NORMr
l = geo.l / geo.r;
% raggio rotore NORMr
x = geo.xr / geo.r;
% lunghezza gap NORMr
g = geo.g / geo.r;
% lunghezza del dente NORMr
lt = geo.lt / geo.r;
% spessore del giogo NORMr
ly = 1 - (x + g) - lt;
% Area denti NORMrpi
vdenti=2*kt*b*x*lt;  
% corrente totale ^2 NORM cfr.eq.(3.3.19) pag. 3.28
% area cave norm. alla circonferenza esterna (NORMrpi)
io2pu = lt*(2*(1-ly)-lt)-vdenti; 
iopu = sqrt(io2pu);
% lunghezza delle testate NORMr
% formula empirica -> prof. Tassoni
ltestata=(2*lt+kracc*(0.5*pi*(1-ly+x)/p)); 
% lunghezza NORMl di pacco attivo e testate
lprimor=(l+ltestata)./l; 
% resistività del rame alla temperatura (TetaCu)
rocu = 17.8*(234.5+tempcu)/(234.5+20)*1e-9; 
% cfr.eq.(3.3.8) pag.3.26
k=sqrt(kcu*kj/(rocu*lprimor));
% Coeff. di NORM delle correnti
CoeffCorrente=pi*(geo.r/1000)^1.5*k/(3*n_spire);  
% MODULO DI CORRENTE
io=CoeffCorrente.*iopu;



