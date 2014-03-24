function assign_block_prop_rot(BLKLABELS,fem,group)
BLKLABELSrot=BLKLABELS.rotore;

% pulisce le selezioni precedenti
mi_clearselected

% Assegna aria alle barriere di flux:
for kk=1:size(BLKLABELSrot.names.BarName,1)
%     keyboard
    mi_addblocklabel(BLKLABELSrot.xy(kk,1),BLKLABELSrot.xy(kk,2));
    mi_selectlabel(BLKLABELSrot.xy(kk,1),BLKLABELSrot.xy(kk,2));
    mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(kk,3)}, 0, fem.res,'None', 0, group, 0);
    mi_clearselected;
end

index=kk+1;
% Assegna ferro di rotore
    mi_addblocklabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
    mi_selectlabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
    mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(index,3)}, 0, fem.res,'None', 0, group, 0);
    mi_clearselected;
    
index=index+1;
% Assegna materiale albero
    mi_addblocklabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
    mi_selectlabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
    mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(index,3)}, 0, fem.res,'None', 0, group, 0);
    mi_clearselected;

end