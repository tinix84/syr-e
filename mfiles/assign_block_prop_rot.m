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
        
        for kk = 1:size(BLKLABELSrot.xy,1)
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
                if geo.BarFillFac == 2
                    angle = atan2(BLKLABELSrot.xy(kk,2),BLKLABELSrot.xy(kk,1))*180/pi+180;
                end
                mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(kk,3)}, 0, fem.res, 'None',angle, group, 0);
                if countdx==geo.dx
                    countdx=1;
                    countPM=countPM+1;
                else
                    countdx=countdx+1;
                end
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
            tmp = [mat.LayerMag.Br mat.LayerMag.Br];
            Br = repmat(tmp,1,geo.ps);
            % Br=[Br 0*Br];
        end
%         for kk=1:length(Br)
%             if (Br(kk)==0)
%                 % air
%                 mi_addblocklabel(BLKLABELSrot.xy(kk,1),BLKLABELSrot.xy(kk,2));
%                 mi_selectlabel(BLKLABELSrot.xy(kk,1),BLKLABELSrot.xy(kk,2));
%                 mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(kk,3)}, 0, fem.res,'None', 0, group, 0);
%                 mi_setblockprop('Air', 0, fem.res,'None', 0, group, 0);
%                 mi_clearselected;
%             else
%                 % Plasto-Magnet with assigned Br
%                 Hc=1/(4e-7*pi)*Br(kk);        % Propriet?da assegnare al magnete
%                 magdir=atan2(BLKLABELSrot.xy(kk,7),BLKLABELSrot.xy(kk,6))*180/pi;
%                 mi_addmaterial([mat.LayerMag.MatName '_' num2str(kk)], 1, 1, Hc);
%                 mi_addblocklabel(BLKLABELSrot.xy(kk,1),BLKLABELSrot.xy(kk,2));
%                 mi_selectlabel(BLKLABELSrot.xy(kk,1),BLKLABELSrot.xy(kk,2));
%                 mi_setblockprop([mat.LayerMag.MatName '_' num2str(kk)], 0, fem.res,'None', magdir, group, 0);
%                 mi_clearselected;
%             end
%         end
%         
%         % rotor iron
%         index=kk+1;
%         mi_addblocklabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
%         mi_selectlabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
%         mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(index,3)}, 0, fem.res,'None', 0, 22, 0);
%         mi_clearselected;
%         
%         % shaft
%         index=index+1;
%         mi_addblocklabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
%         mi_selectlabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
%         mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(index,3)}, 0, fem.res,'None', 0, group, 0);
%         mi_clearselected;
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
                        Hc=1/(4e-7*pi)*Br(kk); 
                        magdir=atan2(BLKLABELSrot.xy(ii,7),BLKLABELSrot.xy(ii,6))*180/pi;
                        mi_addmaterial([mat.LayerMag.MatName '_' num2str(kk)], 1, 1, Hc);
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
                    
        
end
