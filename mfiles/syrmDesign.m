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
%   one machine will be saved and visualized in syre

clc

mu0 = 4e-7*pi;                                                              % air permeability
% DATA not (yet) included into syre GUI
Bs = 2.4;                                                                   % saturation flux density in the ribs [T]

% load('mot_01.mat');

Bfe = dataSet.Bfe;                                                          % steel loading (yoke flux density [T])
kt = dataSet.kt;                                                            % kt = wt/wt_unsat

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

flag_dx=1;      % flag_dx=0 --> dx=0
% flag_dx=1 --> constant rotor carrier width
% flag_dx=2 --> rotor carrier width proportional to sine integral
% flag_dx=3 --> rotor carrier width proportional to flux of d-axis staircase (kt needed)



switch dataSet.TypeOfRotor
    case 'SPM'
        % design domain according to x and lm
        m = 31; n = 21;                                                     % m x n grid of evaluated machines
        lm_g = linspace(dataSet.bRange(1),dataSet.bRange(2),m);
        x = linspace(dataSet.xRange(1),dataSet.xRange(2),n);                % rotor/stator split factor
        ap = geo.phi/180;                                                   % PM/pole ratio

        % parametric analysis: design domain (x,b)
        [xx,lm_gp] = meshgrid(x,lm_g);
        
        r = R * xx;
        rocu = 17.8*(234.5 + tempcuest)/(234.5+20)*1e-9;                    % resistivity of copper [Ohm m]
        ssp = (r+g) * pi/(3*p*q);                                           % stator slot pitch (x)
        sso = ssp * acs;                                                    % stator slot opening (x)
        Ar = r*(sqrt(2)-1);                                                 % shaft radius (x) [mm]
        [kw, ~] = calcKwTh0(avv,6*p*q,p);                                   % winding factor calculation
        
        % Bg calculation
        %         lm = lm_gp*g;
        %         xPMco = r;
        %         xPMregular = r-lm + beta*lm; yPMregular = 0;
        %         [xPMo,yPMo] = rot_point(xPMregular,yPMregular,phi/2*pi/180);        % PM edge point
        %         xArccenter = (xPMco + xPMo - (yPMo.^2./(xPMco-xPMo)))/2;            % find arc center location
        %         Rc = r - xArccenter;
        %
        %         csi = linspace(-phi*pi/180/2,phi*pi/180/2,300);
        %         air = zeros(size(csi,1),ceil((180*size(csi,2)/phi-size(csi,2))/2)); % the size of no mag zone relates to Am
        %         for mm = 1:m
        %             for nn = 1:n
        %                 Lm{mm,nn} = (r(mm,nn)-Rc(mm,nn))*cos(csi) + sqrt(Rc(mm,nn)^2-(r(mm,nn)*sin(csi)-Rc(mm,nn)*sin(csi)).^2)-r(mm,nn)+lm(mm,nn);
        %                 G{mm,nn} = lm(mm,nn) +geo.g - Lm{mm,nn};
        %                 kc{mm,nn} = ssp(mm,nn)./(ssp(mm,nn)-2/pi*G{mm,nn}.*(sso(mm,nn)./G{mm,nn}.*atan(sso(mm,nn)./(2*G{mm,nn}))-log(1+(sso(mm,nn)./(2*G{mm,nn})).^2)));
        %                 Bg{mm,nn} = Lm{mm,nn}./G{mm,nn}./(Lm{mm,nn}./G{mm,nn}+kc{mm,nn}*mur)*Br;
        %                 temp{mm,nn} = [air,Bg{mm,nn},air];
        %                 Bg_avg(mm,nn) = mean(Bg{mm,nn});
        %                 Bg_pole{mm,nn} = [temp{mm,nn},-temp{mm,nn}];
        %                 L = length(Bg_pole{mm,nn});
        %                 Y = fft(Bg_pole{mm,nn});
        %                 P2 = abs(Y/L);
        %                 Bg1(mm,nn) = 2*P2(2);                                       % get Bg1 from fft of Bg
        %             end
        %         end
        
        if Br ==0
            h = errordlg('Please use a real magnet material and define Br in Other Options tab');
            uiwait(h);
            return
        end
                
        [Bt_max,Bg_avg,Bg1] = evalBgapSyrmDesign_SPM(geo,mat,lm_g,x);
 
        wt = 2*pi*r.*Bt_max/(6*p*q)/Bfe;
        ly = pi*r.*Bg_avg*ap/(2*p)/Bfe;                                     % Bianchi 'Theory and design of fractional-slot pm machines'(7.1)
        lt = R - r -g - ly;                                                 % slot length (x,lm) [mm]
        lt(lt<2*ttd)=NaN;
        
        % d axis
        Fmd = 2*r.*Bg1'*l*Ns*kw/p*1e-6;
        
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
                [xLT2,yLT2]=intersezione_retta_circonferenza(0,0,x6,mm2,q2);
                mm1 = (yLT2-y6)./(xLT2-x6);                                 % slope of slot bottom line
                mm2 = (y3-yLT2)./(x3-xLT2);                                 % slope of slot side line
                angle1 = atan(abs((mm2-mm1)./(1+mm1.*mm2)));                % angle between two lines (minor than 90)
                area_corner(ii,jj) = RaccordoFC^2 * (1./tan(angle1/2)-(pi-angle1)/2);           % redundant area at bottom slot area
                
                xArea{ii,jj} = [x2 x2 x3 xLT2 x6 x2];
                yArea{ii,jj} = [ 0 y2 y3 yLT2  0  0];
                area_half_slot(ii,jj) = polyarea(xArea{ii,jj},yArea{ii,jj}) - area_corner(ii,jj);
            end
        end
        Aslots = 2 * area_half_slot *6*p*q;
        Aslots(Aslots<0)=NaN;
        
        kj = Loss/(2*pi*R*l)*1e6;                                           % specific loss (x,lm) [W/m2]
        K = sqrt(kcu*kj/rocu*l./(l+lend));                                  % factor K [] (x,lm)
        i0 = pi/(3*Ns)*(R/1000)^1.5*K.*sqrt(Aslots/(pi*R^2));               % rated current i0 [A] pk
        loadpu = dataSet.CurrLoPP;                                          % current load in p.u. of i0
        
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
        Ld = Lmd + Lslot + Ltip;
        Lq = Ld;
        PF = Fmd./sqrt(Fmd.^2 + (0.001*Lq.* iq).^2);
        
        figure(1)
        [c, h] = contour(x,lm_g,T); clabel(c,h,'FontSize',12); grid on, hold on
        [c, h] = contour(x,lm_g,PF,0.7:0.01:1); clabel(c,h,'FontSize',12); grid on,
        xlabel('x - rotor / stator split'), ylabel('lm/g - magnet / airgap split');
        legend('Torque [Nm]','PF');
        
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
            geo.r = geo.x*R;                                                % rotor radius [mm]
            geo.wt = interp2(x,lm_g,wt,geo.x,geo.lm_g);                     % tooth width [mm]
            geo.wt = round(geo.wt*100)/100;
            geo.lt=interp2(x,lm_g,lt,geo.x,geo.lm_g);                       % slot length [mm]
            geo.lt=round(geo.lt*100)/100;
            geo.Ar=interp2(x,lm_g,Ar,geo.x,geo.lm_g);                       % shaft radius [mm]
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
        m = 31; n = 21;                                                     % m x n grid of evaluated machines
        b = linspace(dataSet.bRange(1),dataSet.bRange(2),m);                % iron/copper split factor
        x = linspace(dataSet.xRange(1),dataSet.xRange(2),n);                % rotor/stator split factor
        
        %         [~, ~, geo, per, ~] = data0(dataSet);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % parametric analysis: design domain (x,b)
        [xx,bb] = meshgrid(x,b);
        
        kt = 1;
        
        [xGap,yGap,kf1,kfm] = evalBgapSyrmDesign(q,kt);                     % airgap induction shape, first harmonics factor and mean value factor
        
        r = R*xx;
        rocu = 17.8*(234.5 + tempcuest)/(234.5+20)*1e-9;                    % resistivity of copper [Ohm m]
        ssp = r * pi/(3*p*q);                                               % stator slot pitch (x,b)
        sso = ssp * acs;                                                    % stator slot opening (x,b)
        kc = ssp./(ssp-2/pi*g*(sso/g.*atan(sso/(2*g))-log(1+(sso/(2*g)).^2)));   % Carter coefficient (x,b)
        ly = pi/2*R/p*xx.*bb*kfm;                                           % yoke or back iron (x,b) [mm]
        %ly = R/p*xx.*bb;                                                   % yoke [mm], do not depend on kt
        wt = 2*pi*R/(6*p*q)*xx.*bb.*kt;                                     % tooth width (x,b) [mm]
        cos_x0 = cos(pi/2/p);
        sin_x0 = sin(pi/2/p);
        Ar = geo.R*xx*(1/cos_x0-sqrt(((1-cos_x0^2)/cos_x0)^2+sin_x0^2));    % (max) shaft radius (x,b) [mm]
        
        lt = R*(1-xx)-g-ly;
        lt(lt<2*ttd) = NaN;
        
        % total insulation la
        % 0) la(x,b) = r - Ar - ly (total radial space r - Ar: too much insulation)
        % la0 = R * xx - Ar - ly;
        % 1) la(x,b) = x0 - rbeta - ly (reduced radial space x0-rbeta-Ar: too much iron)
        % la1 = x0 - rbeta - Ar - ly;
        x0 = r /cos(pi/2/p);                                                % center of barriers circles
        geo.dalpha = geo.dalpha_pu*(90/p);                                  % [mec degrees]
        beta_temp = atand(r*sind(geo.dalpha(1))./(x0 - r * cosd(geo.dalpha(1))));
        rbeta = (x0 - r*cosd(geo.dalpha(1)))./(cosd(beta_temp));            % radius of barrier circle
        % 2) 1st carrier takes 1-cos(p*alpha1) p.u. flux
        la = (x0 - rbeta - Ar - ly*cosd(p*geo.dalpha(1)))*nlay/(nlay-0.5);
        %%
        % rotor design + slot evaluation + Lfqpu + Lcqpu
        alpha = cumsum(geo.dalpha);                                         % alpha in syre coordinates (mech deg, zero is axis Q)
        [df,da] = staircaseAnyAlpha(alpha*p*pi/180);                        % rotor staircase: set of rotor slots positions
        f = cumsum((df));                                                   % stator MMF staircase
        sumDf2r = sum(df.^2);
        Lcqpu = 1-4/pi*sum(f.^2.*da);
        
        alpha = cumsum(da);                                                 % alpha defined as in Vagati tutorial (0 = d axis)
        geo.alpha = cumsum(geo.dalpha);                                     % alpha in syre coordinates (mech deg, zero is axis Q)
        
        
        % initializing matrix
        % q-axis
        Lfqpu = zeros(m,n);
        % rotor
        hc = cell(m,n);
        sk = cell(m,n);
        hf = cell(m,n);
        dx = cell(m,n);
        hc_pu = cell(m,n);
        Bx0 = cell(m,n);
        delta = zeros(m,n);
        lr = zeros(m,n);
        % slot
        xArea = cell(m,n);
        yArea = cell(m,n);
        area_half_slot = zeros(m,n);
        area_corner = zeros(m,n);
        d1 = zeros(m,n);
        c0 = zeros(m,n);
        c1 = zeros(m,n);
        c2 = zeros(m,n);
        
        dfQ=fliplr(df); % df using syre conventions (alpha=0 is the q-axis)
        %% calculate slot area (regualr region subtract redundant region of fillet radius)
        for rr=1:m
            for cc=1:n
                % slot area evaluation
                alpha_slot=2*pi/(ns*p);          % angolo di mezzo passo cava
                RSI=r(rr,cc)+g;                  % r traferro statore
                
                % r eq for middle of slot computation
                mr=tan(alpha_slot/2);
                % tooth side computation
                % design like a line parallel to r, case of trapezoidal slot r2: y=m2x+q2
                % explicit form
                mm2=mr; q2= -wt(rr,cc)/2*sqrt(1+mr^2);
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
                x6=RSI+lt(rr,cc);
                y6=0;
                % LT2 position at the tooth
                [xLT2,yLT2]=intersezione_retta_circonferenza(0,0,x6,mm2,q2);
                mm1 = (yLT2-y6)./(xLT2-x6);       % slope of slot bottom line
                mm2 = (y3-yLT2)./(x3-xLT2);       % slope of slot side line
                angle1 = atan(abs((mm2-mm1)./(1+mm1.*mm2)));        % angle between two lines (minor than 90)
                area_corner(rr,cc) = RaccordoFC^2 * (1./tan(angle1/2)-(pi-angle1)/2);           % redundant area at bottom slot area
                
                xArea{rr,cc} = [x2 x2 x3 xLT2 x6 x2];
                yArea{rr,cc} = [ 0 y2 y3 yLT2  0  0];
                area_half_slot(rr,cc) = polyarea(xArea{rr,cc},yArea{rr,cc}) - area_corner(rr,cc);
                
                % parameters for the evaluation of Ls
                d1(rr,cc)=x3-x2;
                c0(rr,cc)=2*y2;
                c1(rr,cc)=2*y3;
                c2(rr,cc)=2*((xLT2-x6)^2+yLT2^2)^0.5;
                
                % rotor design
                % sk
                beta = calc_apertura_cerchio(geo.alpha*pi/180,geo.R*xx(rr,cc),x0(rr,cc));
                rbeta = (x0(rr,cc) - R*xx(rr,cc) * cosd(geo.alpha))./(cos(beta));
                [xpont,ypont] = calc_intersezione_cerchi(R*xx(rr,cc)-geo.pont0, rbeta, x0(rr,cc));
                sk{rr,cc}=rbeta.*beta;
                lr(rr,cc)=mean(sk{rr,cc}(end-1:end));                       % length of the inner flux carrier (for saturation factor)
                % hc (flux barrier design)
                switch flag_pb
                    case 0 % hc = cost
                        hc{rr,cc}=la(rr,cc)/geo.nlay.*ones(1,geo.nlay);
                        %disp('flux barrier design: hc = cost')
                    case 1 % pb = cost
                        hc{rr,cc}=la(rr,cc)/sum(sk{rr,cc}).*sk{rr,cc};
                        %disp('flux barrier design: pbk = hc/sk = cost')
                    case 2 % min Lfq
                        hc{rr,cc}=la(rr,cc)/sum(dfQ.*sk{rr,cc}).*(dfQ.*sk{rr,cc});
                        %disp('flux barrier desig: hc/(df*sk) = cost')
                end
                rpont_x0=sqrt(ypont.^2+(x0(rr,cc)-xpont).^2);
                Bx0{rr,cc}=x0(rr,cc)-(rpont_x0);
                hc_min=(R*xx(rr,cc)-Ar(rr,cc)-(geo.R-geo.R*xx(rr,cc)-geo.g-lt(rr,cc)))/geo.nlay/4;
                hfe_min=2*geo.pont0;
                delta(rr,cc)=(0.5*hc{rr,cc}(1)+sum(hc{rr,cc}(2:end-1))+0.5*hc{rr,cc}(end)-hc_min*(geo.nlay-1))/(Bx0{rr,cc}(1)-Bx0{rr,cc}(end)-hfe_min*(geo.nlay-1)-hc_min*(geo.nlay-1));
                hc_pu{rr,cc}=hc{rr,cc}*(delta(rr,cc)*geo.nlay)/sum(hc{rr,cc});
                % dx (flux carrier design)
                alphad=[0 90-fliplr(geo.alpha)*geo.p 90];                   % 0<=alphad<=90 [° elt]
                r_all=geo.R*xx(rr,cc);
                for ii=1:geo.nlay
                    r_all=[r_all Bx0{rr,cc}(ii)+hc{rr,cc}(ii)/2 Bx0{rr,cc}(ii)-hc{rr,cc}(ii)/2];
                end
                r_all=[r_all Ar(rr,cc)];
                hf0=r_all(1:2:end)-r_all(2:2:end);
                switch flag_dx
                    case 0
                        hf{rr,cc}=hf0;
                    case 1 % constant iron
                        hf_cost=[ones(1,geo.nlay-1) 0.5]/(geo.nlay-0.5);
                        hf{rr,cc}=hf_cost*sum(hf0(2:end));
                        %disp('flux carrier design: Fe = cost')
                    case 2 % iron proportional to first harmonic flux
                        level=(cosd(alphad(1:end-1))+cosd(alphad(2:end)))/2.*(alphad(2:end)-alphad(1:end-1));
                        level_pu=level/sum(level);
                        hfTemp=fliplr(level_pu)*sum(hf0);
                        hf{rr,cc}=hfTemp(2:end);
                        %disp('flux carrier design: Fe proportional to first harmonic flux')
                    case 3 % iron proportional to flux
                        level_pu=evalSatStairCase(xGap,yGap,alphad);
                        hf{rr,cc}=fliplr(level_pu)*sum(hf0);
                        %disp('flux carrier design: Fe proportional to flux')
                end
                for ii=geo.nlay:-1:1
                    if ii==geo.nlay
                        B1tmp=Ar(rr,cc)+hf{rr,cc}(ii);
                    else
                        B1tmp=B2k(ii+1)+hf{rr,cc}(ii);
                    end
                    dxtmp=1-2/hc{rr,cc}(ii)*(Bx0{rr,cc}(ii)-B1tmp);
                    B1ktmp=Bx0{rr,cc}(ii)-hc{rr,cc}(ii)/2+dxtmp*hc{rr,cc}(ii)/2;
                    B2ktmp=Bx0{rr,cc}(ii)+hc{rr,cc}(ii)/2+dxtmp*hc{rr,cc}(ii)/2;
                    if B1ktmp>(Bx0{rr,cc}(ii)-geo.pont0)
                        B1new=Bx0{rr,cc}(ii)-geo.pont0;
                        dxtmp=1-2/hc{rr,cc}(ii)*(Bx0{rr,cc}(ii)-B1new);
                    elseif B2ktmp<(Bx0{rr,cc}(ii)+geo.pont0)
                        B2new=Bx0{rr,cc}(ii)+geo.pont0;
                        dxtmp=2/hc{rr,cc}(ii)*(B2new-Bx0{rr,cc}(ii))-1;
                    end
                    dx{rr,cc}(ii)=dxtmp;
                    B1k(ii)=Bx0{rr,cc}(ii)-hc{rr,cc}(ii)/2+dx{rr,cc}(ii)*hc{rr,cc}(ii)/2;
                    B2k(ii)=Bx0{rr,cc}(ii)+hc{rr,cc}(ii)/2+dx{rr,cc}(ii)*hc{rr,cc}(ii)/2;
                end
                if flag_dx==0
                    dx{rr,cc}=zeros(1,geo.nlay);
                    %disp('flux carrier design: dx=0')
                    hf{rr,cc}=hf0;
                end
                
                % q-axis flow-through inductance [pu]
                Lfqpu(rr,cc) = 4/pi*p*g*kc(rr,cc)/(xx(rr,cc)*R)*(sum((dfQ).^2.*(sk{rr,cc}./hc{rr,cc})));
            end
        end
        
        %%
        % d axis
        if flag_kw
            [kw, ~] = calcKwTh0(avv,6*p*q,p);
        else
            kw = pi/2/sqrt(3);
        end
        
        Ht = interp1(mat.Stator.BH(:,1),mat.Stator.BH(:,2),Bfe./kt);
        Hy = interp1(mat.Stator.BH(:,1),mat.Stator.BH(:,2),Bfe);
        %ks = 1+mu0*(Hy*2*pi/(6*p*q)*(R-ly/2)+Ht*lt+Hy*mean(sk{rr,cc}(end-1:end)))./(bb.*Bfe.*kf1.*kc.*g);
        ks = 1+mu0*pi/2*(Hy*pi/(6*p*q)*(R-ly/2) + Ht*lt + Hy*lr)./(bb*Bfe*kf1.*kc*g);
        
        Fmd = 2*(R*1e-3)*(l*1e-3)*kw*Ns*Bfe.*kf1/p.*xx.*bb;                 % flux linkage [Vs]
        id = pi*Bfe*kc*(g*1e-3)*p.*ks/(mu0*3*kw*Ns).*bb;                    % id [A]
        
        Lmd = Fmd./id;                                                      % d magnetization inductance Lmd [Vs]
        
        % q axis
        kdq = 1 - Lcqpu - Lfqpu;                                            % anisotrophy factor
        
        % stator design
        if q<1
            lend = 0.5*(wt+pi*(r+lt/2)*sin(pi/(6*p*q)));
        else
            lend = 2*lt+(0.5*pi*(R-ly+r)/p);                                % end turn length (x,b) [mm]
        end
        
        %% calculate slot area and rated current
        Aslots = 2 * area_half_slot *6*p*q;
        
        kj = Loss/(2*pi*R*l)*1e6;                                           % specific loss (x,b) [W/m2]
        K = sqrt(kcu*kj/rocu*l./(l+lend));                                  % factor K [] (x,b)
        
        i0 = pi/(3*Ns)*(R/1000)^1.5*K.*sqrt(Aslots/(pi*R^2));               % rated current i0 [A] pk
        i0=real(i0);
        loadpu = dataSet.CurrLoPP;                                          % current load in p.u. of i0
        
        id(id>loadpu*i0)=i0(id>loadpu*i0);
        gamma=acos(id./(loadpu*i0));
        iq=i0.*sin(gamma);                                                  % q-axis current [A] pk
        
        %iq = sqrt((loadpu*i0).^2 - id.^2); iq = real(iq);                  % q-axis current [A] pk
        Am = Aslots.* id./(loadpu*i0);                                      % slots area dedicated to id [mm2]
        
        Nbob  = Ns/p/(q)/2;                                                 % conductors in slot per layer
        J = 2*Nbob * loadpu*i0 ./ (Aslots/(q*6*p)*kcu);                     % current density in copper [A/mm2] pk
        A = 2*Nbob * loadpu*i0 ./ (r*2*pi/(q*6*p));                         % linear current density [A/mm] pk
        
        % tangential ribs effect
        
        Lrib = 4/pi*kw*Ns*2*(pont0*1e-3)*(l*1e-3)*Bs./iq;
        Lrqpu = Lrib./Lmd;
        
        kdq = kdq-Lrqpu;
        
        
        T = 3/2*p*(kdq.*Fmd.*iq);
        
        % stator leakage inductance Ls
        [dfs] = staircaseRegular(6*q);                                      % stator staircase
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
        
        Lspu = 4/pi*p*kc*g*sumDf2s.*ps./r;                                  % Ls/Lmd: p.u. leakage inductance
        % power factor
        Ld = Lmd.*(1 +Lspu);
        Lq = (Lcqpu + Lfqpu + Lspu + Lrqpu).*Lmd;
        
        csi = Ld./Lq;
        gamma = atand(iq./id);                                              % current phase angle [deg]
        delta = atand((Lq.*iq)./(Ld.*id));                                  % flux linkage phase angle [deg]
        
        PF = sind(gamma-delta);                                             % PF @ gamma (same gamma as torque)
        PFmax = (Ld-Lq)./(Ld+Lq);                                           % PF @ max PF gamma
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
            geo.r = geo.x*R;                                                % rotor radius [mm]
            %     bt = geo.b;                                               % Bgap/Bfe,tooth (tooth p.u. size)
            geo.wt = interp2(x,b,wt,geo.x,geo.b);                           % tooth width
            geo.wt = round(geo.wt*100)/100;
            geo.lt=interp2(x,b,lt,geo.x,geo.b);                             % slot length
            geo.lt=round(geo.lt*100)/100;
            geo.x0=R * geo.x /cos(pi/2/p);
            geo.Ar=interp2(x,b,Ar,geo.x,geo.b);                             % shaft radius [mm]
            geo.Ar=round(geo.Ar*100)/100;
            geo.la=interp2(x,b,la,geo.x,geo.b);                             % total insulation
            
            % hc evaluation - flux barriers design
            geo.alpha=cumsum(geo.dalpha);
            switch flag_pb
                case 0 % hc = cost
                    disp('flux barrier design: hc = cost')
                case 1 % pb = cost
                    disp('flux barrier design: pbk = hc/sk = cost')
                case 2 % min Lfq
                    disp('flux barrier desig: hc/(df*sk) = cost')
            end
            
            switch flag_dx
                case 0 % dx=0
                    disp('flux carrier design: dx=0')
                case 1 % constant iron
                    disp('flux carrier design: Fe = cost')
                case 2 % iron proportional to first harmonic flux
                    disp('flux carrier design: Fe proportional to first harmonic flux')
                case 3 % iron proportional to flux
                    disp('flux carrier design: Fe proportional to flux')
            end
            
            geo.hc_pu = zeros(1,geo.nlay);
            geo.dx = zeros(1,geo.nlay);
            hcTmp=zeros(m,n);
            dxTmp=zeros(m,n);
            for ii=1:geo.nlay
                for mm=1:m
                    for nn=1:n
                        hcTmp(mm,nn)=hc_pu{mm,nn}(ii);
                        dxTmp(mm,nn)=dx{mm,nn}(ii);
                    end
                end
                geo.hc_pu(ii)=interp2(xx,bb,hcTmp,geo.x,geo.b);
                geo.dx(ii)=interp2(xx,bb,dxTmp,geo.x,geo.b);
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
            dataSet.HCpu=round(geo.hc_pu*100)/100;
            dataSet.DepthOfBarrier=round(geo.dx*100)/100;
            
            button = questdlg('pick up another machine?','SELECT','Yes','No','Yes');
            
            if isequal(button,'No')
                buttonS = questdlg('save the last machine?','SELECT','Yes','No','Yes');
                figure(hfig), hold off
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
