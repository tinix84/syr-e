function FemmProblem=assign_block_prop_rotX(FemmProblem,BLKLABELS,fem,group)
BLKLABELSrot=BLKLABELS.rotore;

% Assegna aria alle barriere di flux:
if not(isempty(BLKLABELSrot.names.BarName))
    for kk=1:size(BLKLABELSrot.names.BarName,1)
        %mi_addblocklabel(BLKLABELSrot.xy(kk,1),BLKLABELSrot.xy(kk,2));
        %mi_selectlabel(BLKLABELSrot.xy(kk,1),BLKLABELSrot.xy(kk,2));
        %mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(kk,3)}, 0, fem.res,'None', 0, group, 0);
        %mi_clearselected;
        FemmProblem = addblocklabel_mfemm(FemmProblem, BLKLABELSrot.xy(kk,1),BLKLABELSrot.xy(kk,2), ...
            'BlockType', BLKLABELS.materials{BLKLABELSrot.xy(kk,3)}, ...
            'MaxArea', fem.res,...
            'InGroup', group,...
            'Turns',1);
    end
else
    kk=0;
end
index=kk+1;
% Assegna ferro di rotore
%     mi_addblocklabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
%     mi_selectlabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
%     mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(index,3)}, 0, fem.res,'None', 0, group, 0);
%     mi_clearselected;
    FemmProblem = addblocklabel_mfemm(FemmProblem, BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2), ...
                                  'BlockType', BLKLABELS.materials{BLKLABELSrot.xy(index,3)}, ...
                                  'MaxArea', fem.res,...
                                  'InGroup', group,...
                                  'Turns',0);
    
index=index+1;
% Assegna materiale albero
%     mi_addblocklabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
%     mi_selectlabel(BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2));
%     mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(index,3)}, 0, fem.res,'None', 0, group, 0);
%     mi_clearselected;
    FemmProblem = addblocklabel_mfemm(FemmProblem, BLKLABELSrot.xy(index,1),BLKLABELSrot.xy(index,2), ...
                                  'BlockType', BLKLABELS.materials{BLKLABELSrot.xy(index,3)}, ...
                                  'MaxArea', fem.res,...
                                  'InGroup', group,...
                                  'Turns',0);

end