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

%% Calcolo kavv e la fase della corrente th0
function [kavv, Asse_Fase1] = calcKwTh0(Strati)

% INPUT
% Strati     % Distribuzione degli avvolgimenti (una periodicita' completa anche se la macchina puo' essere simulata sfruttando l'antiperiodicita')

% OUTPUT
% kavv       % Fattore di avvolgimento
% Asse_Fase1 % Posizione dell'asse asse della fase 1 misurato in gradi elettrici rispetto all'asse x (si ipotizza di disegnare la prima cava di statore a cavallo dell'asse x)

% Star of slots

[nLayers,nSlots] = size(Strati);

Passo_elt = 2*pi/nSlots;
Posiz_cave = zeros(1,nSlots);
for ee = 2:nSlots
    Posiz_cave(ee) = Posiz_cave(ee-1) + Passo_elt;
    while Posiz_cave(ee)>2*pi
        Posiz_cave(ee)=Posiz_cave(ee)-2*pi;
    end
end

% Seleziono i vettori della stella degli avvolgimenti appartenenti alla fase 1 e ne eseguo la somma vettoriale
% La fase del vettore risultante identifica la posizione dell'asse della fase 1
num=0; den=0;
for tt=1:nLayers
    for yy=1:nSlots
        if Strati(tt,yy)==1 || Strati(tt,yy)==-1
            den=den+1;
            num = num+Strati(tt,yy)*exp(1i*(Posiz_cave(yy)+0.5*pi));
        end
    end
end

% Fattore di avvolgimento
kavv = abs(num)/den;
% Posizione dell'asse della fase 1 (radianti elettrici) rispetto all'asse x
Asse_Fase1 = (angle(num));
Asse_Fase1 = Asse_Fase1*180/pi;