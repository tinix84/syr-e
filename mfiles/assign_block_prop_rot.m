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
% RotType = geo.RotType;
Q = geo.ns*geo.p;                    % number of slots

switch geo.RotType
    case 'SPM'
        
        countPM=1;
        countdx=1;
        for kk = 1:size(BLKLABELSrot.xy,1);
            mi_addblocklabel(BLKLABELSrot.xy(kk,1),BLKLABELSrot.xy(kk,2));
            mi_selectlabel(BLKLABELSrot.xy(kk,1),BLKLABELSrot.xy(kk,2));
            if (isnan(BLKLABELSrot.xy(kk,6)))
                mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(kk,3)}, 0, fem.res,'None', 0, group, 0);
            else
                if geo.PMdir=='p'
                    angle = 180/geo.p/2+180/geo.p*(countPM-1)+180*rem(countPM,2);
                else
                    angle = atan2(BLKLABELSrot.xy(kk,2),BLKLABELSrot.xy(kk,1))*180/pi+180*rem(countPM,2);
                end
                
                mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(kk,3)}, 0, fem.res, 'None',angle, group, 0);
                if countdx==geo.dx
                    countdx=1;
                    countPM=countPM+1;
                else
                    countdx=countdx+1;
                end
%                 ii = real(fix((kk-1)/geo.dx));
%                 if ((6*geo.t/Q)>1)
%                     if mod(ii,2)==0
%                         mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(kk,3)}, 0, fem.res,'None', 180+(2*ii+1)*180/geo.t/geo.ps, group, 0);
%                     else
%                         mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(kk,3)}, 0, fem.res,'None', (2*ii+1)*180/geo.t/geo.ps, group, 0);
%                     end
%                 else
%                     if mod(ii,2)==0
%                         mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(kk,3)}, 0, fem.res,'None', 90/geo.t/geo.ps+ii*180/geo.ps+180, group, 0);
%                     else
%                         mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(kk,3)}, 0, fem.res,'None', 90/geo.t/geo.ps+ii*180/geo.ps, group, 0);
%                     end
%                 end
            end
            mi_clearselected;
        end
        % rotor iron
        index=kk-1;
        mi_addblocklabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
        mi_selectlabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
        mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(index,3)}, 0, fem.res,'None', 0, 22, 0);
        mi_clearselected;
    otherwise
        % Assegna aria alle barriere di flux:
        %     for kk=1:length(BLKLABELSrot.BarName)
        if ((geo.BarFillFac~=0)&&strcmp('Circular',geo.RotType))
            Br=[Br 0*Br];
        end
        for kk=1:length(Br)
            if (Br(kk)==0)
                % air
                mi_addblocklabel(BLKLABELSrot.xy(kk,1),BLKLABELSrot.xy(kk,2));
                mi_selectlabel(BLKLABELSrot.xy(kk,1),BLKLABELSrot.xy(kk,2));
                mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(kk,3)}, 0, fem.res,'None', 0, group, 0);
                mi_setblockprop('Air', 0, fem.res,'None', 0, group, 0);
                mi_clearselected;
            else
                % Plasto-Magnet with assigned Br
                Hc=1/(4e-7*pi)*Br(kk);        % Proprietà da assegnare al magnete
                magdir=atan2(BLKLABELSrot.xy(kk,7),BLKLABELSrot.xy(kk,6))*180/pi;
                mi_addmaterial([mat.LayerMag.MatName num2str(kk)], 1, 1, Hc);
                mi_addblocklabel(BLKLABELSrot.xy(kk,1),BLKLABELSrot.xy(kk,2));
                mi_selectlabel(BLKLABELSrot.xy(kk,1),BLKLABELSrot.xy(kk,2));
                mi_setblockprop([mat.LayerMag.MatName num2str(kk)], 0, fem.res,'None', magdir, group, 0);
                mi_clearselected;
            end
        end
        
        % rotor iron
        index=kk+1;
        mi_addblocklabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
        mi_selectlabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
        mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(index,3)}, 0, fem.res,'None', 0, 22, 0);
        mi_clearselected;
        
        % shaft
        index=index+1;
        mi_addblocklabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
        mi_selectlabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
        mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(index,3)}, 0, fem.res,'None', 0, group, 0);
        mi_clearselected;
        
end
