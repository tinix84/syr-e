% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.


% Input: dalpha [deg], hc_pu [p.u.]
% Output: alpha [deg], hc [mm]

function geo = calc_alpha_hc_delta_x0_5(geo)

xr = geo.xr;            % Raggio del rotore al traferro
p = geo.p;              % Paia poli
nlay = geo.nlay;        % N° layers
R = geo.r;              % Raggio ext
g = geo.g;              % Traferro
lt = geo.lt;            % Lunghezza denti
pont0 = geo.pont0;      % Ponticelli al traferro


dalpha = geo.dalpha;
hc_pu = geo.hc_pu;
hfe_min = geo.hfe_min;        % min tickness of each steel flux guide

if strcmp(geo.RotType,'ISeg_HS')<1
    pont0=pont0*ones(1,geo.nlay);
    hfe_min=hfe_min*ones(1,geo.nlay);
end

% Eval alpha
alpha = integral_fx(1:length(dalpha),dalpha);
geo.alpha = alpha;

% x0 is the coordinate of the center of the circular barriers
x0 = xr/cos(pi/2/p);
geo.x0 = x0;
% max allowed shaft radius
Ar = x0 - xr * tan(pi/2/p);
geo.ArMaxAdmis = Ar;  
% rotor space available radialwise (air + steel)
htot = xr - Ar;
ly = R - xr - g - lt;       % stator yoke
lyr = 1.0 * ly;             % lower limit of the total steel tickness
la = xr - Ar -lyr;          % upper limit of the total air insulation
% 2014/02/24 MG determination of the minimum thickness of air
% length...
    hc_half_min = la/nlay/8;      % occhio che nn deve essere troppo piccolo se no le barriere verranno sempre eccessivamente piccole, ma?! :-|
%     hc_half_min=0.5*pont0;
 
beta = 180/pi * calc_apertura_cerchio(pi/180*alpha,xr,x0);
r = (x0 - xr * cos(alpha*pi/180))./(cos(beta*pi/180));
% Per il momento la taratura è a mano:

delta=(1/(nlay)*sum(hc_pu));
hfeqMax=r(end)-r(1)-(nlay-1)*2*hc_half_min;
hfeqMin=(nlay-1)*hfe_min;
% hfeq=hfeqMax-(hfeqMax-hfeqMin)/0.8*(delta-0.2);
hfeq=hfeqMax-(hfeqMax-hfeqMin)*(delta);

la=r(end)-r(1)-hfeq;
coeff_gamma=hc_pu(1)/2+hc_pu(nlay)/2+sum(hc_pu(2:nlay-1));
% la=r(end)-r(1)-(nlay-1)*hfe_min;
% 
% laprimo=la/(nlay-1);

hc = [];
error_type1=1;
conta=1;
% while (sum(error_type1)==1)
%     if conta==5
%         break
%     end
if (nlay==1)
    
    %% max hc according to alpha min
    hc_half_max1 = (alpha*pi/180/(1+alpha*pi/180)*(xr-pont0));
    % (needs division by 2 .. don't know why but it works)
    hc_half_max1 = hc_half_max1 * 2;
    
    %% max hc according to alpha max (27 Jan 2011)
    temp_alpha_hfemin = hfe_min/xr; % rad
    temp_alpha_hc_2 = pi/(2*p) - alpha*pi/180 - temp_alpha_hfemin;
    hc_half_max2 = (temp_alpha_hc_2/(1+temp_alpha_hc_2)*(xr-pont0));
    hc_half_max = min(hc_half_max1,hc_half_max2);
%     hc_pu(1) = 1;
    hc(1) = hc_pu(1) * hc_half_max * 2;
    if hc(1)<2*hc_half_min
        hc(1)=2*hc_half_min;
    end
    
    if hc(1)>2*hc_half_max
        hc(1)=2*hc_half_max;
    end
    
else
    % layer one is the outermost, layer nlay is the innermost
    for jj = nlay:-1:1          
        
        if (jj == nlay)
            
%             hc_half_max = min((x0-r(end)-Ar-hfe_min),0.5*la/(0.5+(sum(hc_pu(2:nlay-1)))/(hc_pu(nlay))+hc_pu(1)/hc_pu(nlay)/2));
            hc_half_max = 0.5*la/(0.5+(sum(hc_pu(2:nlay-1)))/(hc_pu(nlay))+hc_pu(1)/hc_pu(nlay)/2);
            hc_nlay_temp_half=0.5*la/(0.5+(sum(hc_pu(2:nlay-1)))/(hc_pu(nlay))+hc_pu(1)/hc_pu(nlay)/2);
            
                if strcmp(geo.RotType,'ISeg_HS')>0
                   hc(jj) = 2*hc_half_max(jj);
                else
                    hc(jj) = 2*hc_half_max(1); 
                end
                
                if hc(jj)<2*hc_half_min
                hc(jj)=2*hc_half_min;
                end
            
        else
            
            if jj == 1
                
%                 hc_half_max = min((r(jj+1)-r(jj)-0.5*hc(jj+1)-hfe_min),(xr-pont0-x0+r(jj)));
                hc_half_max = min((hc_pu(jj)/hc_pu(nlay))*hc_nlay_temp_half,(xr-pont0-x0+r(jj)));
%                 hc_half_max = min(0.5*hc_pu(jj)*laprimo,(xr-pont0-x0+r(jj)));
                if strcmp(geo.RotType,'ISeg_HS')
                   hc(jj) = 2*hc_half_max(jj);
                else
                    hc(jj) = 2*hc_half_max(1); 
                end
                
                if hc(jj)<2*hc_half_min && hc_half_min<=(xr-pont0(jj)-x0+r(jj))
                    hc(jj)=2*hc_half_min;
                end
            else
%                 hc_half_max =min([(hc_pu(jj)/hc_pu(nlay))*hc_nlay_temp_half,(r(jj+1)-r(jj)-hfe_min),(r(jj)-r(jj-1)-hfe_min)]);
                hc_half_max =(hc_pu(jj)/hc_pu(nlay))*hc_nlay_temp_half;
                if strcmp(geo.RotType,'ISeg_HS')
                   hc(jj) = 2*hc_half_max(jj);
                else
                    hc(jj) = 2*hc_half_max(1); 
                end
                
                if hc(jj)<2*hc_half_min
                    hc(jj)=2*hc_half_min;
                end
                
            end
            
        end
        
    end
end

% 2014/02/25 MG Determinazione dei punti di barriera sull'asse q:
hcIni=hc;
hc=abs(hcIni);
Bx0=x0-r;
B1k=Bx0-hc./2;
B2k=Bx0+hc./2;

for k=1:nlay-1
    %% #4 vincolo 1-n° layer overlap);
    if (B2k(k+1)>=B1k(k))   % questa condizione vale invece per tutte le barriere
        Dint=B2k(k+1)-B1k(k);
        B2p=B1k(k+1)+(hc(k)+hc(k+1)-Dint-hfe_min(k))/(1+hc_pu(k)/hc_pu(k+1));
        B1p=B2k(k)-(hc(k)+hc(k+1)-Dint-hfe_min(k))/(1+hc_pu(k+1)/hc_pu(k));
        error_type1=1;
        disp('#4 conflict between two barrier')
%         keyboard
        %% #5 vincolo 1-n° intersezione arie, spessore lato barriera<pont0/2 --> equa ripartizione aria ferro');
            % condizione vale nel caso in cui muovendosi non c'è più spazio per l'aria condizione critica, 
            % scelte casuali erronee, ma si privilegia la fattibilità della
            % macchina per proseguire nell'ottimizzazione...
        if ((Bx0(k)<B1p)||(Bx0(k+1)>B2p)|| ((B2p-Bx0(k+1))<hc_half_min) || ((Bx0(k)-B1p)<hc_half_min))  
            B1p=Bx0(k)-(Bx0(k)-Bx0(k+1))/3;
            B2p=Bx0(k+1)+(Bx0(k)-Bx0(k+1))/3;
            disp('#5 low space for iron and air')
        end % end #5
        B2k(k+1)=B2p;
        B1k(k)=B1p;
    else
        disp('');
        error_type1(k)=0;
    end % end #4
end % end for nlay

if (B1k(nlay)<geo.Ar+hfe_min(end))    % questa condizione vale per l'ultima barriera di flux
    B1k(nlay)=geo.Ar+hfe_min(end);
    disp('#6 vincolo n° layer interseca albero')
end

% Condition for high value of beta
[xc_temp,yc_temp]=calc_intersezione_cerchi(xr,r(nlay),x0);
dPointEndBar=calc_distanza_punti([xc_temp,yc_temp],[xr*cos(pi/2/p), xr*sin(pi/2/p)]);
if (dPointEndBar<(Bx0(nlay)-B1k(nlay)))
   B1k(nlay)=Bx0(nlay)-dPointEndBar+hfe_min(nlay)/2;
end
% hc_half1=Bx0-B1k;
% hc_half2=B2k-Bx0;
% hc=2*min(hc_half1,hc_half2);
hc=B2k-B1k;
% la=sum([hc(1)/2,hc(nlay)/2,sum(hc(2:nlay-1))])
% conta=conta+1
% 
% end % end while
%% 2014/02/25 MG Condition fo the first flux barrier if hc is to high and nlay bigger the length of the barrier too smal to be drawn:
if ((xr-pont0(1)-hc(1)/2)<=x0-r(1))
    temp_hc1=2*(xr+r(1)-x0-(1.5)*pont0(1));
    if (temp_hc1>0)
        hc(1)=temp_hc1;
    end
end
geo.hcIni=hcIni;    % For control the intial hc value are saved

%% end seurity control 1th flux barrier
%%
% SALVO hc NELLA STRUTTURA 'geo'.
geo.hc = hc;
%% CALCOLO 'nlay_effettivo', 'delta' E 'nr' E LI SALVO IN GEO.
geo.delta = [alpha(1) diff(alpha) (90/p)-alpha(end)];
geo.B1k=B1k;
geo.B2k=B2k;
geo.nr = ceil(90 ./ geo.delta ) * 2;
