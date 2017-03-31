% Copyright 2016
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

function [dataSet,flagS] = syrmDesign(dataSet)
%   syrmDesign
%   Script for preliminary design of a Synchonous Reluctance Machine (SyRM)

%   The equations follow the the literature of closed-form design of SyRMs.
%   Main reference are this tutorial course notes:
%   Lipo, T. A., et al. "Synchronous reluctance drives tutorial." IEEE-IAS Annual Meeting. 1994
%   Chapter 3, presented by Prof. A. Vagati, is the one to go and look for
%
%   syrmDesign produces a parametric study, function of x and b
%   one (x,b) or (x, lm/g) combination can be selected from the figure
%   one mamachine will be saved and visualized in syre

clc


mu0 = 4e-7*pi;                              % air permeability
% DATA not (yet) included into syre GUI
%Bfe = 1.4;                                  % steel loading (yoke flux density [T])
Bs = 2.4;                                   % saturation flux density in the ribs [T]
%kt = 1.0;                                   % kt = wt/wt_unsat
Bfe = dataSet.Bfe;
kt = dataSet.kt;

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

Br = mat.LayerMag.Br;
mur = mat.LayerMag.mu;
Loss = per.Loss;
tempcuest = per.tempcuest;

%flag for SyR design
flag_kw=1;      % flag_kw=0 --> use Vagati's equations, with kw=pi/(2*sqrt(3))
                % flag_kw=1 --> use the winding factor

flag_pb=1;      % flag_pb=0 --> hc         = costant
                % flag_pb=1 --> hc/sk      = costant (reduce harmonics content)
                % flag_pb=2 --> hc/(df*sk) = costant (reduce Lfq)

flag_dx=3;      % flag_dx=0 --> dx=0
                % flag_dx=1 --> constant rotor carrier width
                % flag_dx=2 --> rotor carrier width proportional to sine integral
                % flag_dx=3 --> rotor carrier width proportional to flux of d-axis staircase (kt needed)



switch dataSet.TypeOfRotor
    case 'SPM'
        % design domain according to x and lm
        m = 31; n = 21;                                                 % m x n grid of evaluated machines
        lm_g = linspace(dataSet.bRange(1),dataSet.bRange(2),m);
        x = linspace(dataSet.xRange(1),dataSet.xRange(2),n);            % rotor/stator split factor
                
        % parametric analysis: design domain (x,b)
        [xx,lm_gp] = meshgrid(x,lm_g);
        
        r = R * xx;
        rocu = 17.8*(234.5 + tempcuest)/(234.5+20)*1e-9;               % resistivity of copper [Ohm m]
        ssp = r * pi/(3*p*q);                          % stator slot pitch (x)
        sso = ssp * acs;                                            % stator slot opening (x)
        kc = ssp./(ssp-2/pi*g*(sso/g.*atan(sso/(2*g))-log(1+(sso/(2*g)).^2)));   % Carter coefficient (x)
        s = r*(sqrt(2)-1);                                       % shaft radius (x) [mm]
        [kw, ~] = calcKwTh0(avv,6*p*q,p);               % winding factor calculation
        % Bg calculation
        
        if Br ==0
            h = errordlg('Please use a real magnet material and define Br in Other Options tab');
            uiwait(h);
            return
        end
                                                              % permeanbility
        Cphi = 1;                                                       % ratio between airgap area and PM area
        Bg = Br * lm_gp*g/mur./(lm_gp*g/mur +Cphi* kc * g); % TAB 8 WEMPEC 2015
        
        %%  Rounded PM 07/12/2016
        %% calculate Bg1
        if geo.dx < 1
            Bg = Bg;                             % rounded PM: Bg1 is equivalent with Bg of longest magnet
        else
            Bg = 4/pi*Bg*sind(phi/2);        % rectangle PM : Bg1 is calculated by Fourier transform
        end
        %%  Hybrid PM 08/12/2016
        %         Br2 = 0.5;
        %         Bg2 = Br2 * lm_gp*g/mur./(lm_gp*g/mur +Cphi* kc * g);
        %         Bg = 4/pi*(Bg*cos(7/24*pi)+Bg2*(cos(pi/12)-cos(7/24*pi)));
        %
        if q < 1
            ly = 0.94*r.*Bg/Bfe/p;                           % Soong 2014 ECCE ty
            wt = 4*sqrt(2)*r.*Bg/(6*p*q)/Bfe;            % tooth width (x,lm) [mm]
        else
            ly = 1.2*r.*Bg/Bfe/p;                            % Soong 2014 ECCE ty
            wt = 6*r.*Bg/(6*p*q)/Bfe;                    % tooth width (x,lm) [mm]
        end
        lt = R - r -g - ly;                            % slot length (x,lm) [mm]
        lt(lt<2*ttd)=NaN;
        
        % d axis
        Fmd = pi*l*Ns/(sqrt(3)*p)* Bg.*r*1e-6;       % flux linkage Fmd = Lmd*id [Vs]
        % stator design
        if q<1
            lend = 2*lt + 0.5*(wt+pi*(r+lt/2)*sin(pi/(6*p*q)));
        else
            lend=(2*lt+(0.5*pi*(R-ly+r)/p));             % end turn length (x,lm) [mm]
        end
        
        %% calculate slot area (regualr region subtract redundant region of fillet radius) 
        for ii=1:m
            for jj=1:n
                alpha_slot=2*pi/(ns*p);          % angolo di mezzo passo cava
                RSI=r(ii,jj)+g;                  % r traferro statore
                
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
                [xLT2,yLT2]=intersezione_retta_circonferenza(0,0,x6,mm2,q2);
                mm1 = (yLT2-y6)./(xLT2-x6);       % slope of slot bottom line
                mm2 = (y3-yLT2)./(x3-xLT2);       % slope of slot side line
                angle1 = atan(abs((mm2-mm1)./(1+mm1.*mm2)));        % angle between two lines (minor than 90)
                area_corner(ii,jj) = RaccordoFC^2 * (1./tan(angle1/2)-(pi-angle1)/2);           % redundant area at bottom slot area
                
                xArea{ii,jj} = [x2 x2 x3 xLT2 x6 x2];
                yArea{ii,jj} = [ 0 y2 y3 yLT2  0  0];
                area_half_slot(ii,jj) = polyarea(xArea{ii,jj},yArea{ii,jj}) - area_corner(ii,jj);
            end
        end
        Aslots = 2 * area_half_slot *6*p*q;
        Aslots(Aslots<0)=NaN;
        
        kj = Loss/(2*pi*R*l)*1e6;                           % specific loss (x,lm) [W/m2]
        K = sqrt(kcu*kj/rocu*l./(l+lend));                  % factor K [] (x,lm)
        
        i0 = pi/(3*Ns)*(R/1000)^1.5*K.*sqrt(Aslots/(pi*R^2));                       % rated current i0 [A] pk
        loadpu = dataSet.CurrLoPP;                                      % current load in p.u. of i0
        
        %% Inductance calculation
        
        % megnetization inductance
        % use the formula 3.110 of Pyrhonen: total Lm
        Lmd = (6/pi*mu0)/p^2*r.*l./(lm_gp*g+g*kc).*(Ns*kw)^2;   %[mH] Juha Pyrhonen 'Design of rotating electrical machines' (3.110)
        
        % slot leakage inductance, dependent on slot shape
        h1 = ttd;
        b1 = (r+g)*pi/(3*p*q)*acs;
        h2 = 0;
        h4 = lt-h1-h2;
        b4 = pi*(r+g+0.5*lt)./(3*p*q)-wt;
        Lslot = 12*(h4./b4/3+h1./b1)/(6*p*q)*mu0*l*(Ns*kw)^2;                     %[mH] Juha Pyrhonen ' Design of rotating electrical machines'
        
        % tip leakage inductance
        if q >1
            t = 1-kracc;
        else
            t = 1-p/Qs;
        end
        Ltip = 5*((g+lm_gp*g/mur)./b1)./(5+4*(g+lm_gp*g/mur)./b1)*3*4/(6*p*q)*mu0*l*(Ns*kw)^2*(1-3/4*t);  % [mH] Juha Pyrhonen ' Design of rotating electrical machines'
        iq = i0*loadpu;                                                        % q-axis current [A] pk
        T = 3/2*p*Fmd.*iq;
        
        %% PF calculation
        Ld = Lmd + Lslot + Ltip;
        Lq = Ld;
        PF = Fmd./sqrt(Fmd.^2 + (0.001*Lq.* iq).^2);
        
        figure(1)
        [c, h] = contour(x,lm_g,T); clabel(c,h,'Color','Red','FontSize',28); grid on, hold on
        [c, h] = contour(x,lm_g,PF,0.7:0.01:1); clabel(c,h,'FontSize',28); grid on,
        xlabel('x - rotor / stator split'), ylabel('lm/g - magnet / airgap split');
        legend('Torque [Nm]','PF');
        %         title('torque and PF tradeoff');
        
        button = questdlg('pick up a machine?','SELECT','Yes','No','Yes');
        while isequal(button,'Yes')
            figure(1)
            [geo.x,geo.lm_g] = ginput(1);
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            disp(['x = ' num2str(geo.x) ';']);
            disp(['lm_g = ' num2str(geo.lm_g) ';']);
            disp(['Torque = ' num2str(interp2(xx,lm_gp,T,geo.x,geo.lm_g)) ' Nm;']);
            disp(['PwrFac = ' num2str(interp2(xx,lm_gp,PF,geo.x,geo.lm_g))]);
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            
            % Export to Syre GUI
            geo.r = geo.x*R;                                            % rotor radius [mm]
            geo.wt = interp2(x,lm_g,wt,geo.x,geo.lm_g);                 % tooth width [mm]
            geo.wt = round(geo.wt*100)/100;
            geo.lt=interp2(x,lm_g,lt,geo.x,geo.lm_g);                   % slot length [mm]
            geo.lt=round(geo.lt*100)/100;
            geo.Ar=interp2(x,lm_g,s,geo.x,geo.lm_g);                    % shaft radius [mm]
            geo.Ar=round(geo.Ar*100)/100;
            
            % adjourn dataSet
            dataSet.AirGapRadius=round(geo.r*100)/100;
            dataSet.ShaftRadius = round(geo.Ar*100)/100;
            dataSet.ToothLength=geo.lt;
            dataSet.ToothWidth=geo.wt;
            dataSet.ThicknessOfPM = round(geo.lm_g*g*100)/100;
            
            button = questdlg('pick up another machine?','SELECT','Yes','No','Yes');
            
            if isequal(button,'No')
                buttonS = questdlg('save the last machine?','SELECT','Yes','No','Yes');
                figure(1), hold off
            end
            
        end
        
        if ~exist('buttonS')
            buttonS='No';
        end
        
        flagS=0;
        
        if isequal(buttonS,'Yes')
            flagS=1;
            newnamestring = ['x' num2str(geo.x,2) 'lm_g' num2str(geo.lm_g,2)];
            newnamestring(newnamestring=='.') = '';
            dataSet.currentfilename = strrep(dataSet.currentfilename,'.mat',[newnamestring '.mat']);
        end
    otherwise
        % design domain according to b and x
        m = 31; n = 21;                                                 % m x n grid of evaluated machines
        b = linspace(dataSet.bRange(1),dataSet.bRange(2),m);            % iron/copper split factor
        x = linspace(dataSet.xRange(1),dataSet.xRange(2),n);            % rotor/stator split factor
        
%         [~, ~, geo, per, ~] = data0(dataSet);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % parametric analysis: design domain (x,b)
        [xx,bb] = meshgrid(x,b);
        
        [xGap,yGap,kf1,kfm] = evalBgapSyrmDesign(q,kt); % airgap induction shape, first harmonics factor and mean value factor
        
        r = R*xx;
        rocu = 17.8*(234.5 + tempcuest)/(234.5+20)*1e-9;               % resistivity of copper [Ohm m]
        ssp = r * pi/(3*p*q);                          % stator slot pitch (x,b)
        sso = ssp * acs;                                            % stator slot opening (x,b)
        kc = ssp./(ssp-2/pi*g*(sso/g.*atan(sso/(2*g))-log(1+(sso/(2*g)).^2)));   % Carter coefficient (x,b)
        ly = pi/2*R/p*xx.*bb*kfm;                                     % yoke or back iron (x,b) [mm]
        %ly = R/p*xx.*bb;                                        % yoke [mm], do not depend on kt
        wt = 2*pi*R/(6*p*q)*xx.*bb.*kt;               % tooth width (x,b) [mm]
        cos_x0 = cos(pi/2/p);
        sin_x0 = sin(pi/2/p);
        s = geo.R*xx*(1/cos_x0-sqrt(((1-cos_x0^2)/cos_x0)^2+sin_x0^2)); % (max) shaft radius (x,b) [mm]
        
        lt = R*(1-xx)-g-ly;
        lt(lt<2*ttd) = NaN;
        
        % total insulation la
        % 0) la(x,b) = r - s - ly (total radial space r - s: too much insulation)
        % la0 = R * xx - s - ly;
        % 1) la(x,b) = x0 - rbeta - ly (reduced radial space x0-rbeta-s: too much iron)
        % la1 = x0 - rbeta - s - ly;
        x0 = r /cos(pi/2/p);                               % center of barriers circles
        geo.dalpha = geo.dalpha_pu*(90/p);                          % [mec degrees]
        beta_temp = atand(r*sind(geo.dalpha(1))./(x0 - r * cosd(geo.dalpha(1))));
        rbeta = (x0 - r*cosd(geo.dalpha(1)))./(cosd(beta_temp)); % radius of barrier circle
        % 2) 1st carrier takes 1-cos(p*alpha1) p.u. flux
        la = (x0 - rbeta - s - ly*cosd(p*geo.dalpha(1)))*nlay/(nlay-0.5);
        
        % d axis
        if flag_kw
            [kw, ~] = calcKwTh0(avv,6*p*q,p);
        else
            kw=pi/2/sqrt(3);
        end
        
        theta=linspace(0,pi/2,200);
        Hy = 2/pi*trapz(theta,interp1(mat.Stator.BH(:,1),mat.Stator.BH(:,2),Bfe*cos(theta)));
        Ht = interp1(mat.Stator.BH(:,1),mat.Stator.BH(:,2),Bfe);
        clear theta
        ks = 1 + 2*pi/p*mu0*Hy./(kfm*bb*Bfe).*(2*R-ly)./(g*kc) + mu0*Ht./(kfm*bb*Bfe).*lt./(g*kc); % full yoke
        
        
        Fmd = 2*(R*1e-3)*(l*1e-3)*kw*Ns*Bfe.*kf1/p.*xx.*bb;             % flux linkage [Vs] (only 1st harmonic of Bg, mean value)
        id = pi*Bfe*kc*(g*1e-3)*p.*ks/(mu0*3*kw*Ns).*bb.*kfm;           % id [A] (mean value of Bg, all harmonics, but similar to 2/pi*kf1)
        
        Lmd = Fmd./id;                                                  % d magnetization inductance Lmd [Vs]
        
        % q axis
        alpha = cumsum(geo.dalpha);                                     % alpha in syre coordinates (mech deg, zero is axis Q)
        [df,da] = staircaseAnyAlpha(alpha*p*pi/180);                % rotor staircase: set of rotor slots positions
        f = cumsum((df));                                               % stator MMF staircase
        sumDf2r = sum(df.^2);
        Lcqpu = 1-4/pi*sum(f.^2.*da);
        
        Lfqpu = zeros(size(Lcqpu));
        alpha = cumsum(da);                                             % alpha defined as in Vagati tutorial (0 = d axis)
        for j = 1:size(xx,1)
            for k = 1:size(xx,2)
                beta_temp = atan(r(j,k) * sin(alpha) ./ (x0(j,k) - r(j,k) * cos(alpha)));
                rbeta = (x0(j,k) - r(j,k) * cos(alpha))./(cos(beta_temp));
                sk = rbeta .* beta_temp;
                if flag_pb==0
                    hc=la(j,k)/geo.nlay.*ones(1,geo.nlay);
                elseif flag_pb==1
                    hc = la(j,k)/sum(sk)*sk;
                else
                    hc=geo.la/sum(df.*sk).*(df.*sk);
                end
                %Lfqpu(j,k) = 4/pi*p*g*kc(j,k)/la(j,k)*(sum(df.*sqrt(sk/xx(j,k)/R))).^2;
                Lfqpu(j,k) = 4/pi*p*g*kc(j,k)/(xx(j,k)*R)*(sum((df).^2.*(sk./hc)));
            end
        end
        
        kdq = 1 - Lcqpu - Lfqpu;                                        % anisotrophy factor
        
        % stator design
        if q<1
            lend = 2*lt + 0.5*(wt+pi*(r+lt/2)*sin(pi/(6*p*q)));
        else
            lend=(2*lt+(0.5*pi*(R-ly+r)/p));             % end turn length (x,b) [mm]
        end
        
        %% calculate slot area (regualr region subtract redundant region of fillet radius) 
        for ii=1:m
            for jj=1:n
                alpha_slot=2*pi/(ns*p);          % angolo di mezzo passo cava
                RSI=r(ii,jj)+g;                  % r traferro statore
                
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
                [xLT2,yLT2]=intersezione_retta_circonferenza(0,0,x6,mm2,q2);
                mm1 = (yLT2-y6)./(xLT2-x6);       % slope of slot bottom line
                mm2 = (y3-yLT2)./(x3-xLT2);       % slope of slot side line
                angle1 = atan(abs((mm2-mm1)./(1+mm1.*mm2)));        % angle between two lines (minor than 90)
                area_corner(ii,jj) = RaccordoFC^2 * (1./tan(angle1/2)-(pi-angle1)/2);           % redundant area at bottom slot area
                
                xArea{ii,jj} = [x2 x2 x3 xLT2 x6 x2];
                yArea{ii,jj} = [ 0 y2 y3 yLT2  0  0];
                area_half_slot(ii,jj) = polyarea(xArea{ii,jj},yArea{ii,jj}) - area_corner(ii,jj);
                
                % parameters for the evaluation of Ls
                d1(ii,jj)=x3-x2;
                c0(ii,jj)=2*y2;
                c1(ii,jj)=2*y3;
                c2(ii,jj)=2*((xLT2-x6)^2+yLT2^2)^0.5;
            end
        end
        Aslots = 2 * area_half_slot *6*p*q;
        
        kj = Loss/(2*pi*R*l)*1e6;                           % specific loss (x,b) [W/m2]
        K = sqrt(kcu*kj/rocu*l./(l+lend));                  % factor K [] (x,b)
        
        i0 = pi/(3*Ns)*(R/1000)^1.5*K.*sqrt(Aslots/(pi*R^2)); % rated current i0 [A] pk
        loadpu = dataSet.CurrLoPP;                                      % current load in p.u. of i0
        
        iq = sqrt((loadpu*i0).^2 - id.^2); iq = real(iq);               % q-axis current [A] pk
        Am = Aslots.* id./(loadpu*i0);                                  % slots area dedicated to id [mm2]
        
        Nbob  = Ns/p/(q)/2;                             % conductors in slot per layer
        J = 2*Nbob * loadpu*i0 ./ (Aslots/(q*6*p)*kcu); % current density in copper [A/mm2] pk
        A = 2*Nbob * loadpu*i0 ./ (r*2*pi/(q*6*p));  % linear current density [A/mm] pk
        
        % tangential ribs effect
        
        % Fr = 8/pi*kw*Ns*(geo.pont0*1e-3)*(geo.l*1e-3)*Bs;           % rotor ribs flux [Vs]
        % Frpu = 4/pi*p*(pont0/R)*Bs/Bfe;                             % rotor ribs flux [pu]
        Frpu = pi/3*p*pont0/R*Bs./(kf1*Bfe);
        Lqr = Frpu./iq;
        
        kdq = kdq-Lqr;
        
        %T = 3/2*p*(kdq.*Fmd.*iq-Fr.*id);                            % torque (b,x)
        
        T = 3/2*p*(kdq.*Fmd.*iq);
        
        % stator leakage inductance Ls
        [dfs] = staircaseRegular(6*q);                              % stator staircase
        f = cumsum(dfs);
        sumDf2s = sum(dfs.^2);
        d0 = ttd;
        %d1 = d0;
        %d1 = d0*sind(tta);
        d2 = lt- d0 - d1;
        %c0 = acs * R * pi/(6*p*q) * xx;
        %c1 = r * pi/(6*p*q).*(1-bb.*kt);
        %c2 = c1 .* (R - ly)./r;
        beta = c1./c2;
        h = (beta.^2-beta.^4/4-log(beta)-0.75)./((1-beta).*(1-beta.^2).^2);
        ps = d0./c0 + d1./c0.*1./(c1./c0-1).*log(c1./c0)+d2./c2.*h;
        
        Lspu = 4/pi*p*kc*g*sumDf2s.*ps./r;               % Ls/Lmd: p.u. leakage inductance
        % power factor
        Ld = Lmd.*(1 +Lspu);
        Lq = (Lcqpu + Lfqpu + Lspu + Lqr).*Lmd;
        
        csi = Ld./Lq;
        gamma = atand(iq./id);                                          % current phase angle [deg]
        delta = atand((Lq.*iq)./(Ld.*id));                              % flux linkage phase angle [deg]
        
        PF = sind(gamma-delta);                                         % PF @ gamma (same gamma as torque)
        PFmax = (Ld-Lq)./(Ld+Lq);                                       % PF @ max PF gamma
        % PF2 = (csi-1)./sqrt(csi.^2.*(sind(gamma)).^-2+(cosd(gamma)).^-2);  % alternative formula
        
        hfig=figure();
        [c, h] = contour(x,b,T,'Color','r'); clabel(c,h); grid on, hold on
        [c, h] = contour(x,b,PF,0.4:0.02:0.96,'Color','b'); clabel(c,h); grid on,
        % [c, h] = contour(x,b,PFmax,0.7:0.02:0.9); clabel(c, h); grid on,
        xlabel('x - rotor / stator split'), ylabel('b - p.u. magnetic loading');
        legend('[Nm]','PF');
        title('torque and PF tradeoff')
        % set(gca,'FontSize',24)
        
        
        button = questdlg('pick up a machine?','SELECT','Yes','No','Yes');
        
        while isequal(button,'Yes')
            
            figure(hfig)
            [geo.x,geo.b] = ginput(1);
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            disp(['x = ' num2str(geo.x) ';']);
            disp(['b = ' num2str(geo.b) ';']);
            disp(['Torque = ' num2str(interp2(xx,bb,T,geo.x,geo.b)) ' Nm;']);
            disp(['PwrFac = ' num2str(interp2(xx,bb,PF,geo.x,geo.b))]);
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            
            
            % Export to Syre GUI
            geo.r = geo.x*R;                                        % rotor radius [mm]
            %     bt = geo.b;                                           % Bgap/Bfe,tooth (tooth p.u. size)
            geo.wt = interp2(x,b,wt,geo.x,geo.b);                       % tooth width
            geo.wt = round(geo.wt*100)/100;
            geo.lt=interp2(x,b,lt,geo.x,geo.b);                         % slot length
            geo.lt=round(geo.lt*100)/100;
            geo.x0=R * geo.x /cos(pi/2/p);
            geo.Ar=interp2(x,b,s,geo.x,geo.b);                          % shaft radius [mm]
            geo.Ar=round(geo.Ar*100)/100;
            geo.la=interp2(x,b,la,geo.x,geo.b);                         % total insulation
            
            % hc evaluation - flux barriers design
            geo.alpha=cumsum(geo.dalpha);
            beta = 180/pi * calc_apertura_cerchio(pi/180*geo.alpha,geo.r,geo.x0);
            rbeta = (geo.x0 - geo.r * cos(geo.alpha*pi/180))./(cos(beta*pi/180));
            [xpont,ypont] = calc_intersezione_cerchi(geo.r-geo.pont0, rbeta, geo.x0);
            thetabeta=2*abs(atan2(ypont,geo.x0-xpont));
            sk=rbeta.*thetabeta;
            [df,~] = staircaseAnyAlpha(geo.alpha*p*pi/180);
            switch flag_pb
                case 0 % hc = cost
                    hc=geo.la/geo.nlay.*ones(1,geo.nlay);
                    disp('flux barrier design: hc = cost')
                case 1 % pb = cost
                    hc=geo.la/sum(sk).*sk;
                    disp('flux barrier design: pbk = hc/sk = cost')
                case 2 % min Lfq
                    %hc=geo.la/sum(df.*sk).*(df.*sk);
                    hc=geo.la/sum(df.*sk).*(df.*sk);
                    disp('flux barrier desig: hc/(df*sk) = cost')
            end
            rpont_x0=sqrt(ypont.^2+(geo.x0-xpont).^2);
            Bx0=geo.x0-(rpont_x0);
            %la_pu=geo.la*geo.nlay/(geo.r-geo.Ar-2*geo.pont0*(geo.nlay-1));
            hc_min=(geo.r-geo.Ar-(geo.R-geo.r-geo.g-geo.lt))/geo.nlay/4;
            hfe_min=2*geo.pont0;
            delta=(0.5*hc(1)+sum(hc(2:end-1))+0.5*hc(end)-hc_min*(geo.nlay-1))/(Bx0(1)-Bx0(end)-hfe_min*(geo.nlay-1)-hc_min*(geo.nlay-1));
            geo.hc_pu=hc*(delta*geo.nlay)/sum(hc);
            
            % evaluation of dx - flux carriers design
            alphad=[0 90-fliplr(geo.alpha)*geo.p 90]; % 0<=alphad<=90 [° elt]
            r_all=[geo.r];
            for ii=1:geo.nlay
                r_all=[r_all Bx0(ii)+hc(ii)/2 Bx0(ii)-hc(ii)/2];
            end
            r_all=[r_all geo.Ar];
            hf0=r_all(1:2:end)-r_all(2:2:end);
            switch flag_dx
                case 1 % constant iron
                    hf_cost=[ones(1,geo.nlay-1) 0.5]/(geo.nlay-0.5);
                    hf=hf_cost*sum(hf0(2:end));
                    disp('flux carrier design: Fe = cost')
                case 2 % iron proportional to first harmonic flux
                    level=(cosd(alphad(1:end-1))+cosd(alphad(2:end)))/2.*(alphad(2:end)-alphad(1:end-1));
                    level_pu=level/sum(level);
                    hf=fliplr(level_pu)*sum(hf0);
                    hf=hf(2:end);
                    disp('flux carrier design: Fe proportional to first harmonic flux')
                case 3 % iron proportional to flux
                    level_pu=evalSatStairCase(xGap,yGap,alphad);
                    hf=fliplr(level_pu)*sum(hf0);
                    disp('flux carrier design: Fe proportional to flux')
            end
            dx=zeros(1,geo.nlay);
            for ii=geo.nlay:-1:1
                if ii==geo.nlay
                    B1tmp=geo.Ar+hf(ii);
                else
                    B1tmp=B2k(ii+1)+hf(ii);
                end
                dxtmp=1-2/hc(ii)*(Bx0(ii)-B1tmp);
                B1ktmp=Bx0(ii)-hc(ii)/2+dxtmp*hc(ii)/2;
                B2ktmp=Bx0(ii)+hc(ii)/2+dxtmp*hc(ii)/2;
                if B1ktmp>(Bx0(ii)-geo.pont0)
                    B1new=Bx0(ii)-geo.pont0;
                    dxtmp=1-2/hc(ii)*(Bx0(ii)-B1new);
                elseif B2ktmp<(Bx0(ii)+geo.pont0)
                    B2new=Bx0(ii)+geo.pont0;
                    dxtmp=2/hc(ii)*(B2new-Bx0(ii))-1;
                end
                dx(ii)=dxtmp;
                B1k(ii)=Bx0(ii)-hc(ii)/2+dx(ii)*hc(ii)/2;
                B2k(ii)=Bx0(ii)+hc(ii)/2+dx(ii)*hc(ii)/2;
            end
            if flag_dx==0
                dx=0;
                disp('flux carrier design: dx=0')
                hf=hf0;
            end
            geo.dx=dx;
            
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
            dataSet.HCpu=round(geo.hc_pu*100)/100;
            dataSet.DepthOfBarrier=round(geo.dx*100)/100;
            
            button = questdlg('pick up another machine?','SELECT','Yes','No','Yes');
            
            if isequal(button,'No')
                buttonS = questdlg('save the last machine?','SELECT','Yes','No','Yes');
                figure(1), hold off
            end
            
        end
        
        if ~exist('buttonS')
            buttonS='No';
        end
        
        flagS=0;
        
        if isequal(buttonS,'Yes')
            
            % save new machine
            flagS=1;
            newnamestring = ['x' num2str(geo.x,2) 'b' num2str(geo.b,2)];
            newnamestring(newnamestring=='.') = '';
            dataSet.currentfilename = strrep(dataSet.currentfilename,'.mat',[newnamestring '.mat']);
        end
end
