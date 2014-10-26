function [hc,dalpha] = Plot_Machine(h,dataSet,flag_plot)
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

%% Plot GUI

axes(h);
cla(h);

[bounds, geo, per] = data0(dataSet);
RQ = dataSet.RQ;
p = geo.p;                      % paia poli
nlay = geo.nlay;                % numero delle barriere
[geo,gamma] = interpretRQ(RQ,geo);
geo.x0 = geo.r/cos(pi/2/geo.p);
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
t = gcd(round(ns*geo.p),geo.p);  % periodicity
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
fem.res = 0;
fem.res_traf = 0;
%% ROTOR
% build the matrices that describe the rotor
ROTmatr;

dalpha = geo.dalpha;            % Angoli dalpha
hc = geo.hc;                    % Altezze in mm

%% STATOR
% builds the matrices that describe the stator
STATmatr;

%% === PLOT ===============================================================
if strcmp(flag_plot,'Y')
    
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
end

end

