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

%% Plot GUI

h = handles.axes5;
axes(h);
cla(h);

dataSet = handles.dataSet;
[bounds, geo, per] = data1(dataSet);

RQ = dataSet.RQ;

p = geo.p;                      % paia poli
nlay = geo.nlay;                % numero delle barriere

if ~isempty(RQ)
    geo.pathname = cd;
    options.iteration = 0;
    options.currentgen = 1;
    options.PopulationSize = 1;
    first_index = 2;
    last_index = first_index + nlay - 1;
    dalpha_pu = RQ(first_index:last_index);
    % if the sum of the pu angles is too large, it is scaled down
    if sum(dalpha_pu) > 1
        dalpha_pu = dalpha_pu/sum(dalpha_pu);
    end
    % dalpha(2) to dalpha(nlay) in degrees
    dalpha_temp = dalpha_pu * (90/p - RQ(1));
    % all dalpha in mec degrees
    geo.dalpha = [RQ(1) dalpha_temp(1:end-1)];
    % SPESSORE DELLE BARRIERE: 'hc_pu'
    first_index = last_index + 1;
    last_index = first_index + nlay - 1;
    size(RQ);
    getComputerName();
    geo.hc_pu = RQ(first_index:last_index);
    if (strcmp(geo.RotType,'Fluid')||strcmp(geo.RotType,'Seg'))
        first_index = last_index + 1;
        last_index = first_index + nlay - 1;
        geo.Dfe=RQ(first_index:last_index);
    end
    gamma = RQ(end);
    geo.x0 = geo.xr/cos(pi/2/geo.p);
end

ns = geo.ns;
p = geo.p;
% rotor offset angle
th_m0 = 0;                          % [mec deg]
% offset angle: for coordinate transformations
geo.th0 = th_m0*p - 0;     % [elt deg]
%%
%% evaluation of the number of poles to simulate and corresponding periodicity
%%
Q = geo.ns*geo.p;
t = gcd(round(ns)*geo.p,geo.p);  % periodicity
if ((6*t/Q)>1)
    ps = 2*p/t;
    Qs = Q/t;
else
    ps = p/t;
    Qs = Q/2/t;
end
%assign to variable geo, number of poles and slots simulated.
geo.Qs = Qs;
geo.ps = ps;
% Boundary tipology
if (rem(geo.ps,2) == 0)
    periodicity = 4;
else
    periodicity = 5;
end

fem.res = 0;
fem.res_traf = 0;

%% ROTOR
% build the matrixes which describe the rotor
ROTmatr;
%% STATOR
% builds the matrixed which describe the stator
STATmatr;
%% === PLOT ===============================================================

Mat = [rotore; CavaMat];

[nrig,ncol] = size(Mat);
for i = 1 : nrig
    if Mat(i,ncol) == 0
        x1 = Mat(i,3);
        y1 = Mat(i,4);
        x2 = Mat(i,1);
        y2 = Mat(i,2);
        plot(h,[x1,x2],[y1,y2],'Linewidth',2,'Color','k');
        hold on
        grid on
        grid minor
        axis equal
    else
        dati = Mat(i,:);
        % centro
        x0 = dati(1); y0 = dati(2);
        if dati(7) > 0
            % punto 1
            x1 = dati(3); y1 = dati(4);
            %  punto 2
            x2 = dati(5); y2 = dati(6);
        else
            % punto 1
            x2 = dati(3); y2 = dati(4);
            %  punto 2
            x1 = dati(5); y1 = dati(6);
        end
        r = sqrt((x0 - x1)^2 + (y0 - y1)^2);
        ang1 = atan2((y1-y0),(x1-x0));
        ang2 = atan2((y2-y0),(x2-x0));
        if ang1 < 0
            ang1 = ang1 + 2*pi;
        end
        if ang2 < 0
            ang2 = ang2 + 2*pi;
        end
        theta = linspace(ang1,ang2,200);
        if (ang2 - ang1) < 0
            ang1 = ang1 - 2*pi;
            theta = linspace(ang1,ang2,200);
        end
        x_n = r*cos(theta(2:end-1)) + x0;
        y_n = r*sin(theta(2:end-1)) + y0;
        x = [x1 x_n x2];
        y = [y1 y_n y2];
        plot(h,x,y,'Linewidth',2,'Color','k');
        hold on
        grid on
        grid minor
        axis equal
    end
end
hold off









