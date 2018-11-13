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

function geo = assign_block_prop_rot(BLKLABELS,geo,mat,fem,group)

BLKLABELSrot=BLKLABELS.rotore;
Br = repmat(mat.LayerMag.Br,1,geo.ps);

% pulisce le selezioni precedenti
mi_clearselected
Q = geo.ns*geo.p;                    % number of slots

% Assegna aria alle barriere di flux:
if ((max(geo.BarFillFac)~=0)&&strcmp('Circular',geo.RotType))
    tmp = [mat.LayerMag.Br mat.LayerMag.Br];
    Br = repmat(tmp,1,geo.ps);
elseif strcmp(geo.RotType,'SPM')
    Br = mat.LayerMag.Br(1)*ones(size(BLKLABELSrot.xy(:,1)));
end

kk=1; % tiene conto di quale magnete sto assegnando
for ii=1:length(BLKLABELSrot.xy(:,1))
    switch BLKLABELSrot.xy(ii,3)
        case 1  % Aria
            mi_addblocklabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
            mi_selectlabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
            mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(ii,3)}, 0, fem.res,'None', 0, group, 0);
            mi_setblockprop('Air', 0, fem.res,'None', 0, group, 0);
            mi_clearselected;
        case 6 % PM
            if isfield(mat.LayerMag,'BH')
                magdir=atan2(BLKLABELSrot.xy(ii,2),BLKLABELSrot.xy(ii,1))*180/pi;
                mi_addblocklabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                mi_selectlabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                mi_setblockprop(mat.LayerMag.MatName, 0, fem.res,'None', magdir, group, 0);
                mi_clearselected;
            else
                Hc=1/(4e-7*pi*mat.LayerMag.mu)*Br(kk);
                magdir=atan2(BLKLABELSrot.xy(ii,7),BLKLABELSrot.xy(ii,6))*180/pi;
                mi_addmaterial([mat.LayerMag.MatName '_' num2str(kk)], mat.LayerMag.mu, mat.LayerMag.mu, Hc);
                mi_addblocklabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                mi_selectlabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                mi_setblockprop([mat.LayerMag.MatName '_' num2str(kk)], 0, fem.res,'None', magdir, group, 0);
                mi_clearselected;
                kk=kk+1;
            end
        case 5 % Ferro rotore
            mi_addblocklabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
            mi_selectlabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
            mi_setblockprop(mat.Rotor.MatName, 0, fem.res,'None', 0, 22, 0);
            mi_clearselected;
        case 7 % shaft
            mi_addblocklabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
            mi_selectlabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
            if isequal(mat.Shaft.MatName,'ShaftAir')
                mi_setblockprop('Air', 0, fem.res,'None', 0, group, 0);
            else
                mi_setblockprop(mat.Shaft.MatName, 0, fem.res,'None', 0, group, 0);
            end
            mi_clearselected;
    end
end

