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
%   Lipo, T. A., et al. "Synchronous reluctance drives tutorial."�IEEE-IAS Annual Meeting. 1994
%   Chapter 3, presented by Prof. A. Vagati, is the one to go and look for
%
%   syrmDesign produces a parametric study, function of x and b
%   one (x,b) or (x, lm/g) combination can be selected from the figure
%   one machine will be saved and visualized in syre

clc,
[~, ~, geo, per, mat] = data0(dataSet);
%% Analytical model
if strcmp(dataSet.TypeOfRotor,'SPM')
    map = syrmDesign_SPM(dataSet);
elseif strcmp(dataSet.TypeOfRotor,'Vtype')
    map = syrmDesign_Vtype(dataSet);
    warning('Work in progress on Vtype geometry')
else
    map = syrmDesign_SyR(dataSet);
end

%% FEAfix
if dataSet.FEAfixN==0
    map.kd=ones(size(map.xx));
    map.kq=ones(size(map.xx));
    map.xRaw=[];
    map.bRaw=[];
else
    [FEAfixOut]=FEAfix(dataSet,geo,map,dataSet.FEAfixN);
    map.kd   = FEAfixOut.kd;
    map.kq   = FEAfixOut.kq;
    map.xRaw = FEAfixOut.xRaw;
    map.bRaw = FEAfixOut.bRaw;
end

map.fd = map.fd.*map.kd;
map.fq = map.fq.*map.kq;
map.T  = 3/2*geo.p*(map.fd.*map.iq-map.fq.*map.id);
map.PF = abs(sin(atan(map.iq./map.id)-atan(map.fq./map.fd)));

%% Output figure
clc
hfig=figure();
figSetting(15,10)
[c, h] = contour(map.xx,map.bb,map.T,'Color','r','LineWidth',1,'DisplayName','$T$ [Nm]');
clabel(c,h);
[c, h] = contour(map.xx,map.bb,map.PF,0.4:0.02:0.96,'Color','b','LineWidth',1,'DisplayName','$cos \varphi$');
clabel(c,h);
if ~isempty(map.xRaw)
    plot(map.xRaw,map.bRaw,'Color',[0 0.5 0],'LineStyle','none','Marker','o','MarkerFaceColor',[0 0.5 0],'DisplayName','FEAfix')
end
plot(map.xx(isnan(map.T)),map.bb(isnan(map.T)),'rx','DisplayName','unfeasible','MarkerSize',8)
xlabel('$x$ - rotor / stator split');
if strcmp(dataSet.TypeOfRotor,'SPM')
    ylabel('$l_m/g$ - p.u. magnet size')
else
    ylabel('$b$ - p.u. magnetic loading');
end
legend('show','Location','NorthEast')
title('torque and PF tradeoff')

set(hfig,'UserData',map);

%% Machine selection
button = questdlg('pick up a machine?','SELECT','Yes','No','Yes');

while isequal(button,'Yes')
    
    figure(hfig)
    [geo.x,geo.b] = ginput(1);
    if strcmp(dataSet.TypeOfRotor,'SPM')
        setup=inputdlg({'x','lm/g'},'(x,lm/g) values',1,{num2str(geo.x,3),num2str(geo.b,3)});
    else
        setup=inputdlg({'x','b'},'(x,b) values',1,{num2str(geo.x,3),num2str(geo.b,3)});
    end
    if ~isempty(setup)
        geo.x=eval(setup{1});
        geo.b=eval(setup{2});
    end
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp(['x = ' num2str(geo.x) ';']);
    if strcmp(dataSet.TypeOfRotor,'SPM')
        disp(['lm_g = ' num2str(geo.b) ';']);
    else
        disp(['b = ' num2str(geo.b) ';']);
    end
    disp(['Torque = ' num2str(interp2(map.xx,map.bb,map.T,geo.x,geo.b)) ' Nm;']);
    disp(['PwrFac = ' num2str(interp2(map.xx,map.bb,map.PF,geo.x,geo.b))]);
    disp(['CurrDens = ' num2str(interp2(map.xx,map.bb,map.J,geo.x,geo.b)) ' A/mm2']);
    disp(['EltLoading = ' num2str(interp2(map.xx,map.bb,map.A,geo.x,geo.b)) ' A/mm']);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    
    % Export to Syre GUI
    geo.r = geo.x*geo.R;                                % rotor radius [mm]
    geo.wt = interp2(map.xx,map.bb,map.wt,geo.x,geo.b); % tooth width
    geo.wt = round(geo.wt*100)/100;
    geo.lt=interp2(map.xx,map.bb,map.lt,geo.x,geo.b);	% slot length
    geo.lt=round(geo.lt*100)/100;
    if strcmp(dataSet.TypeOfRotor,'SPM')
        dataSet.ThicknessOfPM = round(geo.b*geo.g*100)/100;
    elseif strcmp(dataSet.TypeOfRotor,'Vtype')
        warning('Vtype geometry not yet supported')
        hc_pu = interp2(map.xx,map.bb,map.hc_pu,geo.x,geo.b);
        dataSet.HCpu = round(hc_pu*100)/100;
        beta = interp2(map.xx,map.bb,map.beta,geo.x,geo.b);
        dataSet.SlopeBarrier = round(beta*180/pi*100)/100;
        geo.Ar=interp2(map.xx,map.bb,map.Ar,geo.x,geo.b);                             % shaft radius [mm]
        geo.Ar=round(geo.Ar*100)/100;
    else
        geo.x0=geo.R*geo.x /cos(pi/2/geo.p);
        %geo.Ar=interp2(x,b,Ar,geo.x,geo.b);                             % shaft radius [mm]
        %geo.Ar=round(geo.Ar*100)/100;
        geo.la=interp2(map.xx,map.bb,map.la,geo.x,geo.b);	% total insulation
        
        % hc evaluation - flux barriers design
        geo.alpha=cumsum(geo.dalpha);
        switch map.flag_pb
            case 0 % hc = cost
                disp('flux barrier design: hc = cost')
            case 1 % pb = cost
                disp('flux barrier design: pbk = sk/hc = cost')
            case 2 % min Lfq
                disp('flux barrier desig: hc/(df*sk^0.5) = cost')
        end
        
        switch map.flag_dx
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
        [m,n]=size(map.xx);
        hcTmp=zeros(m,n);
        dxTmp=zeros(m,n);
        for ii=1:geo.nlay
            for mm=1:m
                for nn=1:n
                    hcTmp(mm,nn)=map.hc_pu{mm,nn}(ii);
                    dxTmp(mm,nn)=map.dx{mm,nn}(ii);
                end
            end
            geo.hc_pu(ii)=interp2(map.xx,map.bb,hcTmp,geo.x,geo.b);
            geo.dx(ii)=interp2(map.xx,map.bb,dxTmp,geo.x,geo.b);
        end
        
        dataSet.HCpu=round(geo.hc_pu*100)/100;
        dataSet.DepthOfBarrier=round(geo.dx*100)/100;
    end
    
    % current phase angle
    temp_id = interp2(map.xx,map.bb,map.id,geo.x,geo.b);  % id [A]
    temp_iq = interp2(map.xx,map.bb,map.iq,geo.x,geo.b);  % iq [A]
    dataSet.GammaPP=round(atan2(temp_iq,temp_id)*180/pi*100)/100;
    
    % adjourn dataSet
    dataSet.AirGapRadius=round(geo.r*100)/100;
    dataSet.ShaftRadius = round(geo.Ar*100)/100;
    dataSet.ToothLength=geo.lt;
    dataSet.ToothWidth=geo.wt;
    
    button = questdlg('pick up another machine?','SELECT','Yes','No','Yes');
    
    if isequal(button,'No')
        buttonS = questdlg('save the last machine?','SELECT','Yes','No','Yes');
        figure(hfig)
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

figure(hfig)

if(0)
    
    % Current Density
    figure()
    figSetting(15,10)
    [c, h] = contour(map.xx,map.bb,map.J,'Color','r','LineWidth',1,'DisplayName','$J$ [A/mm2]');
    clabel(c,h);
    [c, h] = contour(map.xx,map.bb,map.A,'Color','b','LineWidth',1,'DisplayName','$A$ [A/mm]');
    clabel(c,h);
    xlabel('$x$ - rotor / stator split');
    if strcmp(dataSet.TypeOfRotor,'SPM')
        ylabel('$l_m/g$ - p.u. magnet size')
    else
        ylabel('$b$ - p.u. magnetic loading');
    end
    title('Current Density Map');
    legend('show','Location','NorthEast')
    
    % Shear Stress
    figure()
    figSetting(15,10)
    
    map.sigma = map.T./(2*pi*geo.R*map.xx*geo.l*1e-6);  % sher stress in Nm/m2
    [c, h] = contour(map.xx,map.bb,map.sigma,'Color','b','LineWidth',1,'DisplayName','$\sigma$ [Nm/m2]');
    clabel(c,h);
    xlabel('$x$ - rotor / stator split');
    if strcmp(dataSet.TypeOfRotor,'SPM')
        ylabel('$l_m/g$ - p.u. magnet size')
    else
        ylabel('$b$ - p.u. magnetic loading');
    end
    title('Shear Stress Vs Cu Temperature Map');
    
    map.tempcuest = map.T*0;
    for j = 1:size(map.T,1)
        for k= 1:size(map.T,2)
            
            geo.r  = map.xx(j,k)*geo.R;
            geo.wt = map.wt(j,k);
            geo.lt = map.lt(j,k);
            geo.ly = map.ly(j,k);
            map.tempcuest(j,k) = temp_est_simpleMod(geo,per);

        end
    end
    
    [c, h] = contour(map.xx,map.bb,map.dTempCu+per.temphous,'Color','r','LineWidth',1,'DisplayName','$\theta$ [$^\circ$C]');
    clabel(c,h);
    legend('show','Location','NorthEast')
        
end
% % 
% figure(1);

