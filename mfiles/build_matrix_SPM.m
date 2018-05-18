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

function rotor = build_matrix_SPM(temp,geo)

seg = geo.dx;
Beta = geo.BarFillFac;
NoSeg = floor(seg);

for jj = 1:floor(NoSeg/2)
    xPMso(jj) = temp.xPMso(jj);
    yPMso(jj) = temp.yPMso(jj);
    xPMsi(jj) = temp.xPMsi(jj);
    yPMsi(jj) = temp.yPMsi(jj);
end
xPMco = temp.xPMco;
yPMco = temp.yPMco;
xPMo = temp.xPMo;
yPMo = temp.yPMo;
xPMci = temp.xPMci;
yPMci = temp.yPMci;
xPMi = temp.xPMi;
yPMi = temp.yPMi;

x4 = temp.x4;
y4 = temp.y4;
x5 = temp.x5;
y5 = temp.y5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if geo.BarFillFac == 2
    
    %% hybrid shape
    xFeo = temp.xFeo;
    yFeo = temp.yFeo;
    xFei = temp.xFei;
    yFei = temp.yFei;
    
    rotor = [0 0 xPMco yPMco xPMo yPMo 1;
        0 0 xPMci yPMci xPMi yPMi 1;
        0 0 xPMo yPMo xFeo yFeo 1;
        0 0 xPMi yPMi xFei yFei 1;
        0 0 xFei yFei x4 y4 1;
        0 0 xFeo yFeo x5 y5 1];
    
    rotor = [ rotor;
        xPMo yPMo xPMi yPMi NaN NaN 0;
        xFeo yFeo xFei yFei NaN NaN 0;
        x4 y4 x5 y5 NaN NaN 0];
else
    %% rectangle or rounded shape
    xArccenter = temp.xArccenter;
    yArccenter = temp.yArccenter;
    x6 = temp.x6;
    y6 = temp.y6;
    
    rotor = [0 0 xPMci yPMci xPMi yPMi 1;
        0 0 xPMi yPMi x4 y4 1];
    %         0 0 x6 y6 x5 y5 1];                         % air box
    
    rotor = [rotor;
        xArccenter yArccenter xPMco yPMco xPMo yPMo 1];
    
    rotor = [ rotor;
        xPMo yPMo xPMi yPMi NaN NaN 0];
    %         x4 y4 x5 y5 NaN NaN 0];
    
%         if Beta <1
%             rotor = [rotor;
%                 x6 y6 xPMo yPMo NaN NaN 0];
%         end
    
end

if seg ~=1
    for jj = 1:floor(NoSeg/2)
        rotor = [rotor; xPMso(jj) yPMso(jj) xPMsi(jj) yPMsi(jj) NaN NaN 0];
    end
end
