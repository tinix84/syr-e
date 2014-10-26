% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

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
%% Data
q=geo.q;
p=geo.p;
% qui dentro si usa notazione normalizzata stile Guglielmi - Vagati
% per cui non badare troppo a significato dei simboli

% superficie esterna statore - m2
Asup = 2*pi*geo.R*geo.l*1e-6; 
% dissipazione - W/m2
kj = per.Loss / Asup; 
% temperatura
tempcu = per.tempcu;
% paia poli
p = geo.p;
% induzione normalizzata nel ferro
% b = geo.b;
% coeff. scaling dente
% kt = geo.kt;
% bt = geo.bt;

b  = (geo.R - geo.r - geo.lt)*geo.p/geo.r;  % Bgap/Bfe,yoke (back-iron p.u. size)
bt = geo.wt/pi*(3*geo.q*geo.p)/geo.r;        % Bgap/Bfe,tooth (tooth p.u. size)
    

% spire in serie per fase
n_spire = geo.Ns;
% raccorciamento di passo
kracc = geo.kracc;
% riempimento cava
kcu = geo.kcu;

% lunghezza pacco NORMr
l = geo.l / geo.R;
% raggio rotore NORMr
x = geo.r / geo.R;
% lunghezza gap NORMr
g = geo.g / geo.R;
% lunghezza del dente NORMr
lt = geo.lt / geo.R;
% spessore del giogo NORMr
ly = 1 - (x + g) - lt;
% Area denti NORMrpi
vdenti=2*bt*x*lt;  
% corrente totale ^2 NORM cfr.eq.(3.3.19) pag. 3.28
% area cave norm. alla circonferenza esterna (NORMrpi)
io2pu = lt*(2*(1-ly)-lt)-vdenti; 
iopu = sqrt(io2pu);
if (geo.q<1)
    wt=2*pi*kt*x*b/(6*p*q);
    ltestata=0.5*(wt+pi*(1+lt/2)*sin(pi/(6*p*q)));
else
    % lunghezza delle testate NORMr
    % formula empirica -> prof. Tassoni
    ltestata=(2*lt+kracc*(0.5*pi*(1-ly+x)/p));
    % lunghezza NORMl di pacco attivo e testate
end
lprimor=(l+ltestata)./l;
% resistività del rame alla temperatura (TetaCu)
rocu = 17.8*(234.5+tempcu)/(234.5+20)*1e-9; 
% cfr.eq.(3.3.8) pag.3.26
k=sqrt(kcu*kj/(rocu*lprimor));
% Coeff. di NORM delle correnti
CoeffCorrente=pi*(geo.R/1000)^1.5*k/(3*n_spire);  
% MODULO DI CORRENTE
io=CoeffCorrente.*iopu;



