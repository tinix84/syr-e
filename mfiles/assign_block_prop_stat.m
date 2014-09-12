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

function assign_block_prop_stat(BLKLABELS,geo,fem,group)
BLKLABELSstat=BLKLABELS.statore;

mi_addcircprop('fase1', 0, 1);
mi_addcircprop('fase2', 0, 1);
mi_addcircprop('fase3', 0, 1);
mi_addcircprop('fase1n', 0, 1);
mi_addcircprop('fase2n', 0, 1);
mi_addcircprop('fase3n', 0, 1);

for kk=1:size(BLKLABELSstat.names.air_slot,1)
    mi_addblocklabel(BLKLABELSstat.xy(kk,1),BLKLABELSstat.xy(kk,2));
    mi_selectlabel(BLKLABELSstat.xy(kk,1),BLKLABELSstat.xy(kk,2));
    mi_setblockprop(BLKLABELS.materials{BLKLABELSstat.xy(kk,3)}, 0, fem.res,'None', 0, group, 1);
    mi_clearselected
end

index=kk+1;
%% avvolgimento
Qs=geo.Qs;
avv=geo.avv;

for jj = 1:Qs
    % strato interno
    for kk=1:size(avv,1)
        
        if avv(kk,jj) > 0
            fase = ['fase' num2str(avv(kk,jj))];
        else
            temp = num2str(avv(kk,jj));
            fase = ['fase' temp(2) 'n'];
        end
        mi_addblocklabel(BLKLABELSstat.xy(index,1),BLKLABELSstat.xy(index,2));
        mi_selectlabel(BLKLABELSstat.xy(index,1),BLKLABELSstat.xy(index,2));
        mi_setblockprop(BLKLABELS.materials{BLKLABELSstat.xy(index,3)}, 0, fem.res, fase, 0, group, 1);
        mi_clearselected
        index=index+1;
    end
end

index=index;

for kk=1:size(BLKLABELSstat.names.FeYoke,1)
    mi_addblocklabel(BLKLABELSstat.xy(index,1),BLKLABELSstat.xy(index,2));
    mi_selectlabel(BLKLABELSstat.xy(index,1),BLKLABELSstat.xy(index,2));
    mi_setblockprop(BLKLABELS.materials{BLKLABELSstat.xy(index,3)}, 0, fem.res,'None', 0, group, 1);
    mi_clearselected;
    index=index+1;
end


end