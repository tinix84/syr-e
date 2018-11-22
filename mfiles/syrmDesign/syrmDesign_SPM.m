function [map] = syrmDesign_SPM(dataSet)
% 
% [map] = syrmDesign_SPM(dataSet)
% 
% Preliminary design for SyR machines. References:
% - Chao Lu ECCE 2017 paper

%% inputs
mu0 = 4e-7*pi;                                                              % air permeability

Bfe = dataSet.Bfe;                                                          % steel loading (yoke flux density [T])

[~, ~, geo, per, mat] = data0(dataSet);

R = geo.R;
p = geo.p;
q = geo.q;
acs = geo.acs;
g = geo.g;
avv = geo.avv;
l = geo.l;
kcu = geo.kcu;
Ns = geo.Ns;
Nbob = geo.Nbob;
kracc = geo.kracc;
Qs = geo.Qs;
phi = geo.phi;
ns = geo.ns;
ttd = geo.ttd;
tta = geo.tta;
RaccordoFC = geo.SFR;
nlay = geo.nlay;
pont0 = geo.pont0;
pont = geo.pont0;%+min(geo.pont);

Br = mat.LayerMag.Br;
mur = mat.LayerMag.mu;
Loss = per.Loss;
tempcu = per.tempcu;

%flags for SyRM design
flag_sd = 1;    % flag_sd=1 --> use subdomain model to compute Fmd,wt,ly
                % flag_sd=0 --> use pure analytical equation to compute Fmd,wt,ly

% design domain according to x and lm
m = 31; n = 21;                                                     % m x n grid of evaluated machines
lm_g = linspace(dataSet.bRange(1),dataSet.bRange(2),m);
x = linspace(dataSet.xRange(1),dataSet.xRange(2),n);                % rotor/stator split factor
ap = geo.phi/180;                                                   % PM/pole ratio

% parametric analysis: design domain (x,b)
[xx,lm_gp] = meshgrid(x,lm_g);

beta = geo.BarFillFac;
r = R * xx;
rocu = 17.8*(234.5 + tempcu)/(234.5+20)*1e-9;                    % resistivity of copper [Ohm m]
ssp = (r+g) * pi/(3*p*q);                                           % stator slot pitch (x)
sso = ssp * acs;                                                    % stator slot opening (x)
Ar = r*(sqrt(2)-1);                                                 % shaft radius (x) [mm]
[kw, ~] = calcKwTh0(avv,6*p*q,p);                                   % winding factor calculation

if ~flag_sd
    %Bg calculation
    lm = lm_gp*g;
    xPMco = r;
    xPMregular = r-lm + beta*lm; yPMregular = 0;
    [xPMo,yPMo] = rot_point(xPMregular,yPMregular,phi/2*pi/180);        % PM edge point
    xArccenter = (xPMco + xPMo - (yPMo.^2./(xPMco-xPMo)))/2;            % find arc center location
    Rc = r - xArccenter;

    csi = linspace(-phi*pi/180/2,phi*pi/180/2,300);
    air = zeros(size(csi,1),ceil((180*size(csi,2)/phi-size(csi,2))/2)); % the size of no mag zone relates to Am
    for mm = 1:m
        for nn = 1:n
            Lm{mm,nn} = (r(mm,nn)-Rc(mm,nn))*cos(csi) + sqrt(Rc(mm,nn)^2-(r(mm,nn)*sin(csi)-Rc(mm,nn)*sin(csi)).^2)-r(mm,nn)+lm(mm,nn);
            G{mm,nn} = lm(mm,nn) +geo.g - Lm{mm,nn};
            kc{mm,nn} = ssp(mm,nn)./(ssp(mm,nn)-2/pi*G{mm,nn}.*(sso(mm,nn)./G{mm,nn}.*atan(sso(mm,nn)./(2*G{mm,nn}))-log(1+(sso(mm,nn)./(2*G{mm,nn})).^2)));
            Bg{mm,nn} = Lm{mm,nn}./G{mm,nn}./(Lm{mm,nn}./G{mm,nn}+kc{mm,nn}*mur)*Br;
            temp{mm,nn} = [air,Bg{mm,nn},air];
            Bg_avg(mm,nn) = mean(Bg{mm,nn});
            Bg_pole{mm,nn} = [temp{mm,nn},-temp{mm,nn}];
            L = length(Bg_pole{mm,nn});
            Y = fft(Bg_pole{mm,nn});
            P2 = abs(Y/L);
            Bg1(mm,nn) = 2*P2(2);                                       % get Bg1 from fft of Bg
        end
    end
    %wt = 2*pi*r.*Bt_max/(6*p*q)/Bfe;
    %ly = pi*r.*Bg_avg*ap/(2*p)/Bfe;                                     % Bianchi 'Theory and design of fractional-slot pm machines'(7.1)
    if q < 1            
        ly = pi*r.*Bg_avg*phi/180/(2*p)/Bfe*(1-acs);                    
        wt = 2*pi*r.*Bg_avg/(6*p*q)/Bfe;                                % Hanselman 'Brushless PM machine design' (9.4)
    else
        ly = pi*r.*Bg_avg*phi/180/(2*p)/Bfe;                            % Hanselman 'Brushless PM machine design' (9.7)
        wt = 2*pi*r.*Bg_avg/(6*p*q)/Bfe;                                % Hanselman 'Brushless PM machine design' (9.4)
    end
else
    [wt,ly,Bg1] = evalBgapSyrmDesign_SPM(geo,mat,lm_g,x,Bfe);
end

if Br ==0
    h = errordlg('Please use a real magnet material and define Br in Other Options tab');
    uiwait(h);
    return
end


%  
lt = R - r -g - ly;                                                 % slot length (x,lm) [mm]
lt(lt<2*ttd)=NaN;

% d axis
if flag_sd
    Fmd = 2*r.*Bg1'*l*Ns*kw/p*1e-6;
else
    Fmd = pi*l*Ns/(sqrt(3)*p)* Bg1.*r*1e-6;
end
% stator design
if q<1
    lend = 0.5*(wt+pi*(r+lt/2)*sin(pi/(6*p*q)));
else
    lend = 2*lt+(r+g+lt/2)*2*pi/(p*q);                              % end-turn length [mm]
end

%% calculate slot area (regualr region subtract redundant region of fillet radius)
for ii=1:m
    for jj=1:n
        alpha_slot=2*pi/(ns*p);                                     % angolo di mezzo passo cava
        RSI=r(ii,jj)+g;                                             % r traferro statore

        % r eq for middle of slot computation
        mr=tan(alpha_slot/2);
        % tooth side computation
        % design like a line parallel to r, case of trapezoidal slot r2: y=m2x+q2
        % explicit form
        mm2=mr; q2= -wt(ii,jj)/2*sqrt(1+mr^2);
        [x1t,y1t]=intersezione_retta_circonferenza(0,0,RSI,mm2,q2);

        slot_open_ang=acs*2*pi/(ns*p)/2;
        [x1,y1]=intersezione_retta_circonferenza(0,0,RSI,tan(slot_open_ang),0);
        if (tan(y1/x1)>tan(y1t/x1t))
            x1=x1t;
            y1=y1t;
        end

        x2=x1+ttd;
        y2=y1;

        mtta=tan(pi/2-tta*pi/180);
        qtta=y2-mtta*x2;
        % ytta=mtta*xx+qtta;
        [x3,y3]=intersezione_tra_rette(mtta,-1,qtta,mm2,-1,q2);

        % end of the slot
        x6=RSI+lt(ii,jj);
        y6=0;
        % LT2 position at the tooth
        if geo.parallel_slot
            q2 = y3;
            m2 = 0;
        end
        [xLT2,yLT2]=intersezione_retta_circonferenza(0,0,x6,mm2,q2);
        mm1 = (yLT2-y6)./(xLT2-x6);                                 % slope of slot bottom line
        mm2 = (y3-yLT2)./(x3-xLT2);                                 % slope of slot side line
        angle1 = atan(abs((mm2-mm1)./(1+mm1.*mm2)));                % angle between two lines (minor than 90)
        area_corner(ii,jj) = RaccordoFC^2 * (1./tan(angle1/2)-(pi-angle1)/2);           % redundant area at bottom slot area

        xArea{ii,jj} = [x2 x2 x3 xLT2 x6 x2];
        yArea{ii,jj} = [ 0 y2 y3 yLT2  0  0];
        area_half_slot(ii,jj) = polyarea(xArea{ii,jj},yArea{ii,jj}) - area_corner(ii,jj);
        
        % copper overtemperature
        geo0=geo;
        geo0.r  = geo.R*xx(rr,cc);
        geo0.wt = wt(rr,cc);
        geo0.lt = lt(rr,cc);
        geo0.ly = ly(rr,cc);
        
        dTempCu(rr,cc) = temp_est_simpleMod(geo0,per);
    end
end

dTempCu=dTempCu-per.temphous;

Aslots = 2 * area_half_slot *6*p*q;
Aslots(Aslots<0)=NaN;

kj = Loss/(2*pi*R*l)*1e6;                                           % specific loss (x,lm) [W/m2]
K = sqrt(kcu*kj/rocu*l./(l+lend));                                  % factor K [] (x,lm)
i0 = pi/(3*Ns)*(R/1000)^1.5*K.*sqrt(Aslots/(pi*R^2));               % rated current i0 [A] pk
loadpu = dataSet.CurrLoPP;                                          % current load in p.u. of i0

Nbob  = Ns/p/q/2;                                                   % conductors in slot per layer
J = 2*Nbob * loadpu*i0 ./ (Aslots/(q*6*p)*kcu);                     % current density in copper [A/mm2] pk
A = 2*Nbob * loadpu*i0 ./ (r*2*pi/(q*6*p));                         % linear current density [A/mm] pk

%% Inductance calculation
kc = ssp./(ssp-2/pi*g*(sso/g.*atan(sso/(2*g))-log(1+(sso/(2*g)).^2)));   % Carter coefficient (x)

% megnetization inductance
% use the formula 3.110 of Pyrhonen: total Lm
Lmd = (6/pi*mu0)/p^2*r.*l./(lm_gp*g+g*kc).*(Ns*kw)^2;               %[mH] Juha Pyrhonen 'Design of rotating electrical machines' (3.110)

% slot leakage inductance, dependent on slot shape
h1 = ttd;
b1 = (r+g)*pi/(3*p*q)*acs;
h2 = 0;
h4 = lt-h1-h2;
b4 = pi*(r+g+0.5*lt)./(3*p*q)-wt;
Lslot_self = 12*(h4./b4/3+h1./b1)/(6*p*q)*mu0*l*(Ns*kw)^2;          %[mH] Juha Pyrhonen ' Design of rotating electrical machines' (4.30)
%% mutual inductance is included
Lslot_mutual = 2*(h4./b4/3+h1./b1)/(6*p*q)*mu0*l*(Ns*kw)^2;
Lslot = Lslot_self + Lslot_mutual;                                  %[mH] El-Refaie "Winding Inductances of Fractional Slot Surface-Mounted Permanent Magnet Brushless Machines," (17)

% tip leakage inductance, short pitching is neglected
Ltip = 5*((g+lm_gp*g/mur)./b1)./(5+4*(g+lm_gp*g/mur)./b1)*3*4/(6*p*q)*mu0*l*(Ns*kw)^2;  % [mH] Juha Pyrhonen ' Design of rotating electrical machines' (4.63)
iq = i0*loadpu;                                                     % q-axis current [A] pk
T = 3/2*p*Fmd.*iq;

%% PF calculation
Ld = (Lmd + Lslot + Ltip)*1e-3;
Lq = Ld;
PF = Fmd./sqrt(Fmd.^2 + (Lq.* iq).^2);


%% save all the results in the map structure
map.xx      = xx;
map.bb      = lm_gp;
map.T       = T;
map.PF      = PF;
map.id      = zeros(m,n);
map.iq      = iq;
map.fd      = Fmd;
map.fq      = Lq.*iq;
map.wt      = wt;
map.lt      = lt;
map.Ar      = Ar;
map.geo     = geo;
map.ly      = ly;
map.Lmd     = Lmd;
map.Ltip    = Ltip;
map.Lslot   = Lslot;
map.Aslots  = Aslots;
map.A       = A;
map.J       = J;
map.dtempCu = dtempCu;
map.dataSet = dataSet;



