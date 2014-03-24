%%
% 2014 MG 

function geo = calc_feseable_hc_alpha(geo)

xr = geo.xr;            % Raggio del rotore al traferro
p = geo.p;              % Paia poli
nlay = geo.nlay;        % N° layers
R = geo.r;              % Raggio ext
g = geo.g;              % Traferro
lt = geo.lt;            % Lunghezza denti
pont0 = geo.pont0;      % Ponticelli al traferro

dalpha = geo.dalpha;
hc_pu = geo.hc_pu;

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
% evaluation of the barriers ticknesses hc
hfe_min = 2 * pont0;        % min tickness of each steel flux guide
%%
%% CONTROLLO DI SUCUREZZA PER IL DISEGNO DELLE BARRIERE (NO SOVRAPPOSIZIONI,... RISP. AMPIEZZE MINIME ecc...)
%%
error_code=[];
Dfe_old=Dfe;
for k=1:nlay  
   if  (hc(k)/2*(1-abs(Dfe(k)))<=pont0);
   disp('#1 Dfe riassegnati');
   Dfe1=1-2*pont0/hc(k);
   error_code=1;
      if (Dfe(k)>0)
       Dfe(k)=Dfe1;
      else
          Dfe(k)=-Dfe1;
      end
   end
end

geo.Dfe=Dfe;
B1k=Bx0-hc/2+Dfe.*hc/2; B2k=Bx0+hc/2+Dfe.*hc/2;

if nlay~=1
    
for k=1:nlay-1
    
    hfemin=2*pont0;
    hcmin=1;
    hc_old=hc;
    if ((xr-B2k(k))<1) % questa condizione varrebbe per la 1°barriera
        B2k(k)=xr-1;
        disp('#2 vincolo 1 layer esce dal rotore');
        error_code=[error_code,2];
    end
    if (B1k(end)<geo.Ar+hfemin)    % questa condizione vale per l'ultima barriera di flux
        B1k(end)=geo.Ar+1.0;
        disp('#3 vincolo n° layer interseca albero')   
        error_code=[error_code,3];
    end                

    if (B2k(k+1)>=B1k(k))   % questa condizione vale invece per tutte le barriere
        Dq=B2k(k+1)-B1k(k);
        B2p=B2k(k+1)-(1/2)*(Dq+hfemin);
        B1p=B1k(k)+(1/2)*(Dq+hfemin);
        disp('#4 vincolo 1-n° layer sovrapposizione'); 
        error_code=[error_code,4];
       
        if ((Bx0(k)<B1p)||(Bx0(k+1)>B2p)|| ((B2p-Bx0(k+1))<pont0/2) || ((Bx0(k)-B1p)<pont0/2))  % condizione vale nel caso in cui muovendosi con #2 non c'è più spazio per l'aria
           B1p=Bx0(k)-(Bx0(k)-Bx0(k+1))/3;
           B2p=Bx0(k+1)+(Bx0(k)-Bx0(k+1))/3;
           disp('#5 vincolo 1-n° intersezione arie, spessore lato barriera<pont0/2 --> equa ripartizione aria ferro');
           error_code=[error_code,5];
        end
        B2k(k+1)=B2p;
        B1k(k)=B1p;
        
    elseif((B1k(k)-B2k(k+1))<hfemin)    %questa condizione vale quando non si è riusciti ad assicurare un ferro minimo tra le barriere di flux (non ho trovato nulla di meglio per adesso, pensarci su!!!!)
        dB12=Bx0(k)-Bx0(k+1);
        Dhc12=dB12-hfemin;
        B1p=Bx0(k)-Dhc12/2;
        B2p=Bx0(k+1)+Dhc12/2;
        disp('#6 vincolo 1-n° tra 2 barriere successive no abbastanza ferro --> equa ripartizione aria ferro');        
        error_code=[error_code,6];
        
        if ((Bx0(k)-B1p)<pont0 || (B2p-Bx0(k+1))<pont0)
        B1p=Bx0(k)-(Bx0(k)-Bx0(k+1))/3;
        B2p=Bx0(k+1)+(Bx0(k)-Bx0(k+1))/3;
        disp('#7 vincolo 1-n° (Bx0(k)-B1p)<pont0 (B2p-Bx0(k+1))<pont0 --> equa ripartizione aria ferro');
        error_code=[error_code,7];        
        end
        
        B2k(k+1)=B2p;
        B1k(k)=B1p;

    end
%
%     if ((Bx0(k)-B1k(k))<(pont0/2))
%         B1k(k)=Bx0(k)+pont0/2;
%                 
%     end
    
end

else
    hfemin=2*pont0;
    hcmin=1;
    hc_old=hc;
    if ((xr-B2k(k))<1) % questa condizione varrebbe per la 1°barriera
        B2k(k)=xr-1;
        disp('#2 vincolo 1 layer esce dal rotore');
        error_code=[error_code,2];
    end

end
error_code=[error_code,zeros(1,3*nlay-length(error_code))];

hc=B2k-B1k; geo.hc=hc;

% Instruction for flux barrier translaction in function of range Dfe=[-1:1]:
for k=1:length(Dfe)
if (Dfe(k)==1)||(Dfe(k)==-1)
    Bx0(k)=(B1k(k)+B2k(k))/2;
    [xBx0New(k),yBy0New(k)]=rot_point(Bx0(k),0,pi/2/p);
    Br0New(k)=abs(xBx0New(k)+1j*yBy0New(k));
    CNew(k)=(((Br0New(k)./rshaft).^(2*p)-1)./(Br0New(k)./rshaft).^p);
    tetaRpontNew(k)=180/p-(1/p)*asin((CNew(k)*((xr-geo.pont0)/rshaft)^p)/(((xr-geo.pont0)/rshaft)^(2*p)-1))*180/pi; 
    % Riassegno xpont ed ypont
    [xpontNew,ypontNew]=rot_point((xr-geo.pont0)*cos(tetaRpontNew(k)*pi/180),(xr-geo.pont0)*sin(tetaRpontNew(k)*pi/180),-pi/2/p); 
    xpont(k)=xpontNew;
    ypont(k)=ypontNew;
end
end

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%
%% SECURITY CONTROL END
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%

end
