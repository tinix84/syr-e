% Copyright 2015
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

function dataSet = syrmDesign(dataSet)
%   syrmDesign
%   Script for preliminary design of a Synchonous Reluctance Machine (SyRM)

%   The equations follow the the literature of closed-form design of SyRMs.
%   Main reference are this tutorial course notes:
%   Lipo, T. A., et al. "Synchronous reluctance drives tutorial." IEEE-IAS Annual Meeting. 1994
%   Chapter 3, presented by Prof. A. Vagati, is the one to go and look for
%
%   syrmDesign produces a parametric study, function of x and b
%   one (x,b) combination can be selected from the figure
%   one mamachine will be saved and visualized in syre

clc

mu0 = 4e-7*pi;          % air permeability
% DATA not (yet) included into syre GUI
Bfe = 1.35;             % steel loading (yoke and tooth flux density [T])


% if eq_set == 1 the design process use the equation of
% 'Boazzo, Pellegrino, Vagati - Multipolar SPM Machines for Direct-Drive
% Application: a General Design Approach'
% for the definition of wt,lt, ly, Lmd, Lslot.
eq_set=0;

if eq_set
    dataSet.bRange(1)=4/pi*sin(dataSet.AngleSpanOfPM/180*pi/2)*dataSet.Br/Bfe*1/(1+1.1*4.5^-1);
    dataSet.bRange(2)=4/pi*sin(dataSet.AngleSpanOfPM/180*pi/2)*dataSet.Br/Bfe*1/(1+1.1*6.5^-1);
end

% design domain according to b and x
m = 31; n = 21;             % m x n grid of evaluated machines
b = linspace(dataSet.bRange(1),dataSet.bRange(2),m);    % iron/copper split factor
x = linspace(dataSet.xRange(1),dataSet.xRange(2),n);    % rotor/stator split factor

[~, ~, geo, per, mat] = data0(dataSet);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parametric analysis: design domain (x,b)
[xx,bb] = meshgrid(x,b);

rocu = 17.8*(234.5 + per.tempcu)/(234.5+20)*1e-9;   % resistivity of copper [Ohm m]
ssp = geo.R * xx * pi/(3*geo.p*geo.q);              % stator slot pitch (x,b)
sso = ssp * geo.acs;                                % stator slot opening (x,b)
kc = ssp./(ssp-2/pi*geo.g*(sso/geo.g.*atan(sso/(2*geo.g))-log(1+(sso/(2*geo.g)).^2)));   % Carter coefficient (x,b)
lt = geo.R * (1-(1+b/geo.p)'*x)-geo.g;              % tooth length (x,b) [mm]
ly = pi/2*geo.R/geo.p*b'*x;                         % yoke or back iron (x,b) [mm]
wt = b' * pi * geo.R * x * 2/(6*geo.p*geo.q);       % tooth width (x,b) [mm]
s = geo.R*xx*(sqrt(2)-1);                           % shaft radius (x,b) [mm]


% from Boazzo - Multipolar SPM
if eq_set
    p=geo.p;
    kt=1;
    Br=mat.LayerMag.Br;
    
    a = pi*geo.R*xx./p;             % pole pitch [mm]
    ly = a/pi.*bb;                  % yoke [mm]
    wt = a/3/geo.q*kt.*bb;          % tooth width [mm]
    lt = geo.R-geo.R*xx-geo.g-ly;   % tooth length [mm]
    lm_g = kc.*(Br/Bfe.*4/pi.*sin(geo.phi/180*pi/2)./bb-1).^-1; %magnet length/airgap length
end

[kw, ~] = calcKwTh0(geo.avv,6*geo.p*geo.q,geo.p); % winding factor isn't equal to 1 for the fractional slot stator and some nonconventional windings

if strcmp(geo.RotType,'SPM')
    if geo.q < 1
        ly = geo.R/geo.p*b'* x;                                   % yoke or back iron (x,b) [mm]   Bianchi's tutorial 2007
        wt = b'* geo.R * x * 4*sqrt(2)/(6*geo.p*geo.q);
    else
        ly =  pi/2*geo.R/geo.p*b'* x*(1-geo.acs);                 % Soong 2014 ECCE ty
    end
    lt = geo.R - geo.R * xx -geo.g - ly;                          % tooth length (x,b) [mm]
end

% total insulation la
% 0) la(x,b) = r - s - ly (total radial space r - s: too much insulation)
% la0 = geo.R * xx - s - ly;
% 1) la(x,b) = x0 - rbeta - ly (reduced radial space x0-rbeta-s: too much iron)
% la1 = x0 - rbeta - s - ly;
x0 = geo.R * xx /cos(pi/2/geo.p);   % center of barriers circles
geo.dalpha = geo.dalpha_pu*(90/geo.p);  % [mec degrees]
beta_temp = atand(geo.R*xx*sind(geo.dalpha(1))./(x0 - geo.R*xx * cosd(geo.dalpha(1))));
rbeta = (x0 - geo.R*xx*cosd(geo.dalpha(1)))./(cosd(beta_temp)); % radius of barrier circle
% 2) 1st carrier takes 1-cos(p*alpha1) p.u. flux
la = (x0 - rbeta - s - ly*cosd(geo.p*geo.dalpha(1)))*geo.nlay/(geo.nlay-0.5);

% d axis
Fmd = pi*geo.R*geo.l*geo.Ns*kw*Bfe/(sqrt(3)*geo.p)*b'*x*1e-6;  % flux linkage Fmd = Lmd*id [Vs]
Bg = Bfe * bb;                                              % peak gap flux density [T]
id = 2*geo.p/(sqrt(3)*geo.Ns)*kc*geo.g/mu0.*Bg*0.001;       % id [A]
Lmd = Fmd./id;                                              % d magnetization inductance Lmd [Vs]

% q axis
alpha = cumsum(geo.dalpha);                             % alpha in syre coordinates (mech deg, zero is axis Q)
[df,da] = staircaseAnyAlpha(alpha*geo.p*pi/180);        % rotor staircase: set of rotor slots positions
f = cumsum((df));                                       % stator MMF staircase
sumDf2r = sum(df.^2);
Lcqpu = 1-4/pi*sum(f.^2.*da);

Lfqpu = zeros(size(Lcqpu));
alpha = cumsum(da);             % alpha defined as in Vagati tutorial (0 = d axis)
for j = 1:size(xx,1)
    for k = 1:size(xx,2)
        beta_temp = atan(geo.R*xx(j,k) * sin(alpha) ./ (x0(j,k) - geo.R*xx(j,k) * cos(alpha)));
        rbeta = (x0(j,k) - geo.R*xx(j,k) * cos(alpha))./(cos(beta_temp));
        sk = rbeta .* beta_temp;
        Lfqpu(j,k) = 4/pi*geo.p*geo.g*kc(j,k)/la(j,k)*(sum(df.*sqrt(sk/xx(j,k)/geo.R))).^2;
    end
end

kdq = 1 - Lcqpu - Lfqpu;        % anisotrophy factor

% stator design
if geo.q<1
    lend = 2*lt + 0.5*(wt+pi*(geo.R*xx+lt/2)*sin(pi/(6*geo.p*geo.q)));
else
    lend=(2*lt+(0.5*pi*(geo.R-ly+xx*geo.R)/geo.p));                         % end turn length (x,b) [mm]
end

Aslots = ((geo.R*xx+lt/2-geo.ttd)*2*pi/(2*6*geo.p*geo.q)-wt/2).*lt*2*6*geo.p*geo.q;         % total slot area (x,b) [mm2] estimate the slot area by using the strcuture in 'Readme'
%     Aslots = pi*geo.R^2*(lt/geo.R.*(2-lt/geo.R-2*b'*x*(1+1/geo.p)));
%     Aslots(Aslots<0) = 0;           % eliminate non feasible solutions (Aslots < 0)

kj = per.Loss/(2*pi*geo.R*geo.l)*1e6;                                                   % specific loss (x,b) [W/m2]
K = sqrt(geo.kcu*kj/rocu*geo.l./(geo.l+lend));                                          % factor K [] (x,b)

i0 = pi/(3*geo.Ns)*(geo.R/1000)^1.5*K.*sqrt(Aslots/(pi*geo.R^2));   % rated current i0 [A] pk
loadpu = dataSet.CurrLoPP;                                          % current load in p.u. of i0

% for SPM
if strcmp(geo.RotType,'SPM')
    
    Br = mat.LayerMag.Br;
    if Br ==0
        h = errordlg('Please use a real magnet material and define Br in Other Options tab');
        uiwait(h);
        return
    end

%     kl = 1;                   % leakage ratio
    mur = 1.1;                % permeanbility
    Cphi = 1;                 % ratio between airgap area and PM area
    lm_general = kc*mur*Cphi*geo.g./(2*Cphi./Bg-1);     % WEMPEC Tab 8 2015
    lm = 2/Br *lm_general;
    
    
    
    %% Inductance calculation
    
    % megnetization inductance
    %Lmd = (8/pi)*(geo.Ns/(2*geo.p))^2*mu0*geo.r*2*geo.l./(lm+geo.g*kc);   % [mH] Juha Pyrhonen ' Design of rotating electrical machines' (3.109)
    %Lmd = (8/pi)*(geo.Ns/(2*geo.p))^2*mu0*geo.R*xx*2*geo.l./(lm+geo.g*kc);   % [mH] Juha Pyrhonen ' Design of rotating electrical machines' (3.109)
    % use the formula 3.110 of Pyrhonen: total Lm
    
    [kw, ~] = calcKwTh0(geo.avv,6*geo.p*geo.q,geo.p); % winding factor isn't equal to 1 for the fractional slot stator and some nonconventional windings
    Lmd = (6/pi*mu0)/(geo.p)^2.*(geo.R*xx).*(geo.l)./(lm+geo.g*kc).*(kw*geo.Ns)^2;  %[mH] Juha Pyrhonen 'Design of rotating electrical machines' (3.110)
    
    
    
    % slot leakage inductance
    h1 = geo.ttd;
    b1 = (geo.r+geo.g)*pi/(3*geo.p*geo.q)*geo.acs;
    % h2 = ((geo.r+geo.g+geo.ttd)*pi/(3*geo.p*geo.q)*0.5-0.5*b1-0.5*geo.wt)*tand(geo.tta);
    h2 = 0;
    h4 = lt-h1-h2;
    b4 = pi*(geo.r+geo.g+0.5*lt)./(3*geo.p*geo.q)-geo.wt;
    Lslot = 12*(h4./b4/3+h1/b1)./(6*geo.p*geo.q)*mu0*geo.l*(geo.Ns)^2;         % [mH] Juha Pyrhonen ' Design of rotating electrical machines'
    
    % tip leakage inductance
    if geo.q >1
        t = 1-geo.kracc;
    else
        t = 1-geo.p/geo.Qs;
    end
    Ltip = 5*((geo.g+lm/mur)/b1)./(5+4*(geo.g+lm/mur)/b1)*3*4/(6*geo.p*geo.q)*mu0*geo.l*(geo.Ns)^2*(1-3/4*t);  % [mH] Juha Pyrhonen ' Design of rotating electrical machines'
    
    
    % from Boazzo - Multipolar SPM
    if eq_set
        q=geo.q;
        ns=length(geo.avv(:,1));    % number of layers
        nr=0;                       % cave di raccorciamento
        Q0 = (6*p*geo.q)./gcd(6*p*geo.q,2*p);
        
        Lbase = mu0*(geo.l/1000)/2*(2/pi*kw*geo.Ns)^2;
        if q>=1
            Lm = pi^2/(6*kw^2)*((1-(q-1)^2/q^3)-nr./(4*q)*(10*q-13-nr*(2*q-1))./(2*q^2))*(1./(kc+lm_g)).*(a./geo.g);
            Ls = pi^2/(2*kw^2)*(lt/geo.g)./(1-bb*kt)*(1-3*nr/16/q).*(geo.g./a);
        else
            Lm = 1/ns*pi^2/(12*(q*kw)^2)*1./(kc+lm_g).*(a./geo.g);
            if strcmp(geo.slot_layer_pos,'over_under')
                Ls = pi^2/(2*kw^2)*(lt/geo.g)./(1-bb*kt).*(1-9*(ns-1)/16/Q0).*(geo.g./a);
            else
                Ls = pi^2/(2*kw^2)*(lt/geo.g)./(1-bb*kt).*(1-3*(ns-1)/4/Q0).*(geo.g./a);
            end
        end
        
        Lmd = Lm*Lbase;
        Lslot = Ls*Lbase;
        
    end
    
    kdq = ones(size(Bg));
    id = zeros(size(Bg));
    

end

iq = sqrt((loadpu*i0).^2 - id.^2); iq = real(iq);                   % q-axis current [A] pk
Am = Aslots.* id./(loadpu*i0);                                      % slots area dedicated to id [mm2]

geo.Nbob  = geo.Ns/geo.p/(geo.q)/2;                                 % conductors in slot per layer
J = 2*geo.Nbob * loadpu*i0 ./ (Aslots/(geo.q*6*geo.p)*geo.kcu);     % current density in copper [A/mm2] pk
A = 2*geo.Nbob * loadpu*i0 ./ (xx*geo.R*2*pi/(geo.q*6*geo.p));      % linear current density [A/mm] pk
T = 3/2*geo.p*kdq.*Fmd.*iq;                                         % torque (b,x)

if eq_set
    B = bb*Bfe;
    A_SPM = 3/2*geo.Ns./(a/1000)*kw.*iq;
    T_SPM = 2*pi*(geo.R*xx/1000).^2.*B.*A_SPM*(geo.l/1000);
end

% stator leakage inductance Ls
[dfs] = staircaseRegular(6*geo.q);                                  % stator staircase
f = cumsum(dfs);
sumDf2s = sum(dfs.^2);
d0 = geo.ttd;
d1 = d0;
d2 = lt- d0 - d1;
c0 = geo.acs * geo.R * pi/(6*geo.p*geo.q) * xx;
c1 = geo.R * pi/(6*geo.p*geo.q)*xx.*(1-bb);
c2 = c1 .* (geo.R - ly)./(geo.R*xx);
beta = c1./c2;
h = (beta.^2-beta.^4/4-log(beta)-0.75)./((1-beta).*(1-beta.^2).^2);
ps = d0./c0 + d1./c1.*1./(c1./c0-1).*log(c1./c0)+d2./c2.*h;

Lspu = 4/pi*geo.p*kc*geo.g/geo.R*sumDf2s.*ps./xx;           % Ls/Lmd: p.u. leakage inductance

% % power factor
Ld = Lmd.*(1 +Lspu);
Lq = (Lcqpu + Lfqpu + Lspu).*Lmd;

% for SPM
if strcmp(geo.RotType,'SPM')
    Ld = Lmd + Lslot +Ltip;
    Lq = Ld;
end

csi = Ld./Lq;
gamma = atand(iq./id);              % current phase angle [deg]
delta = atand((Lq.*iq)./(Ld.*id));  % flux linkage phase angle [deg]

% for SPM
if strcmp(geo.RotType,'SPM')
    delta = atand((Lq.*iq*1e-3)./Fmd);  % flux linkage phase angle [deg]
end

PF = sind(gamma-delta);             % PF @ gamma (same gamma as torque)
PFmax = (Ld-Lq)./(Ld+Lq);           % PF @ max PF gamma
% PF = (csi-1)./sqrt(csi.^2.*(sind(gamma)).^-2+(cosd(gamma)).^-2);  % alternative formula

figure(1)
[c, h] = contour(x,b,T); clabel(c, h); grid on, hold on
[c, h] = contour(x,b,PF,0.6:0.02:0.96); clabel(c, h); grid on,
% [c, h] = contour(x,b,PFmax,0.7:0.02:0.9); clabel(c, h); grid on,
xlabel('x - rotor / stator split'), ylabel('b - p.u. magnetic loading');
legend('[Nm]','PF');
title('torque and PF tradeoff')
% set(gca,'FontSize',24)

% if eq_set
%     figure(2)
%     clf
%     [c, h] = contour(x,b,real(T_SPM)); clabel(c, h); grid on, hold on
%     xlabel('x - rotor / stator split')
%     ylabel('b - p.u. magnetic loading')
%     title('Torque from Multipolar SPM')
% end    


button = questdlg('pick up a machine?','SELECT','Yes','No','Yes');
while isequal(button,'Yes')
    
    figure(1)
    [geo.x,geo.b] = ginput(1);
    display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    display(['x = ' num2str(geo.x) ';']);
    display(['b = ' num2str(geo.b) ';']);
    display(['Torque = ' num2str(interp2(xx,bb,T,geo.x,geo.b)) ' Nm;']);
    display(['PwrFac = ' num2str(interp2(xx,bb,PF,geo.x,geo.b)) ' Nm;']);
    display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    
    button = questdlg('pick up a machine?','SELECT','Yes','No','Yes');
end

button = questdlg('save the last machine?','SELECT','Yes','No','Yes');

if isequal(button,'Yes')
    
    % Export to Syre GUI
    geo.r = geo.x*geo.R;                        % rotor radius [mm]
    %     bt = geo.b;                           % Bgap/Bfe,tooth (tooth p.u. size)
    geo.wt = interp2(x,b,wt,geo.x,geo.b);
    geo.wt = round(geo.wt*100)/100;
    geo.lt=interp2(x,b,lt,geo.x,geo.b);
    geo.lt=round(geo.lt*100)/100;
    geo.x0=geo.R * geo.x /cos(pi/2/geo.p);
    geo.Ar=interp2(x,b,s,geo.x,geo.b);          % shaft radius [mm]
    geo.Ar=round(geo.Ar*100)/100;
    % for SPM
    if strcmp(geo.RotType,'SPM')
        geo.lm=interp2(x,b,lm,geo.x,geo.b);    % PM thickness [mm]
        geo.lm=round(geo.lm*100)/100;
    end
    geo.la=interp2(x,b,la,geo.x,geo.b);        % total insulation
    
    % iterative evaluation of barriers radial dimension hc_pu
    lapu = 1;
    if not(strcmp(geo.RotType,'SPM'))
        geo.hc_pu=lapu*ones(size(geo.dalpha));
        [geo] = calcHcCheckGeoControl(geo);
        lb=sum(geo.hc);                        % Summation of barriers length
        EPS=0.03;
        count=1;
        while abs(geo.la-lb)>=0.1;
%             abs(geo.la-lb)
%             disp(num2str(count))
%             count=count+1;
            if geo.la>lb;
                lapu=lapu*(1+EPS);
            end
            if geo.la<lb;
                lapu=lapu*(1-EPS);
            end
            geo.hc_pu=lapu*ones(1,geo.nlay);
            [geo] = calcHcCheckGeoControl(geo);
            lb = sum(geo.hc);
        end
    end
    
    % current phase angle
    temp_id = interp2(x,b,id,geo.x,geo.b);  % id [A]
    temp_iq = interp2(x,b,iq,geo.x,geo.b);  % iq [A]
    dataSet.GammaPP=round(atand(temp_iq/temp_id)*100)/100;
    
    % adjourn dataSet
    dataSet.AirGapRadius=round(geo.r*100)/100;
    dataSet.ShaftRadius = round(geo.Ar*100)/100;
    dataSet.ToothLength=geo.lt;
    dataSet.ToothWidth=geo.wt;
    dataSet.ThicknessOfPM = geo.lm;
    lapu=round(lapu*100)/100;
    dataSet.HCpu=lapu*ones(1,geo.nlay);
    
    % save new machine
    newnamestring = ['x' num2str(geo.x,2) 'b' num2str(geo.b,2)];
    newnamestring(newnamestring=='.') = '';
    dataSet.currentfilename = strrep(dataSet.currentfilename,'.mat',[newnamestring '.mat']);

end



