function FemmProblem=assign_block_prop_statX(FemmProblem,BLKLABELS,geo,fem,group)
BLKLABELSstat=BLKLABELS.statore;

[FemmProblem, circuitind1] = addcircuit_mfemm(FemmProblem, 'fase1', 'CircType', 1);
[FemmProblem, circuitind2] = addcircuit_mfemm(FemmProblem, 'fase2', 'CircType', 1);
[FemmProblem, circuitind3] = addcircuit_mfemm(FemmProblem, 'fase3', 'CircType', 1);
[FemmProblem, circuitind1n] = addcircuit_mfemm(FemmProblem, 'fase1n', 'CircType', 1);
[FemmProblem, circuitind2n] = addcircuit_mfemm(FemmProblem, 'fase2n', 'CircType', 1);
[FemmProblem, circuitind3n] = addcircuit_mfemm(FemmProblem, 'fase3n', 'CircType', 1);

for kk=1:size(BLKLABELSstat.names.air_slot,1)
    FemmProblem=addblocklabel_mfemm(FemmProblem, BLKLABELSstat.xy(kk,1),BLKLABELSstat.xy(kk,2),...
        'BlockType',BLKLABELS.materials{BLKLABELSstat.xy(kk,3)},...
        'MaxArea',fem.res,...
        'InGroup', group,...
        'Turns',1);
    
    %mi_addblocklabel(BLKLABELSstat.xy(kk,1),BLKLABELSstat.xy(kk,2));
    %mi_selectlabel(BLKLABELSstat.xy(kk,1),BLKLABELSstat.xy(kk,2));
    %mi_setblockprop(BLKLABELS.materials{BLKLABELSstat.xy(kk,3)}, 0, fem.res,'None', 0, group, 1);
    %mi_clearselected
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
        FemmProblem=addblocklabel_mfemm(FemmProblem, BLKLABELSstat.xy(index,1),BLKLABELSstat.xy(index,2),...
            'BlockType',BLKLABELS.materials{BLKLABELSstat.xy(index,3)},...
            'MaxArea',fem.res,...
            'InGroup', group,...
            'InCircuit',fase,...
            'Turns',1);
        %         mi_addblocklabel(BLKLABELSstat.xy(index,1),BLKLABELSstat.xy(index,2));
        %         mi_selectlabel(BLKLABELSstat.xy(index,1),BLKLABELSstat.xy(index,2));
        %         mi_setblockprop(BLKLABELS.materials{BLKLABELSstat.xy(index,3)}, 0, fem.res, fase, 0, group, 1);
        %         mi_clearselected
        index=index+1;
    end
end

index=index;

for kk=1:size(BLKLABELSstat.names.FeYoke,1)
    FemmProblem=addblocklabel_mfemm(FemmProblem, BLKLABELSstat.xy(index,1),BLKLABELSstat.xy(index,2),...
        'BlockType',BLKLABELS.materials{BLKLABELSstat.xy(index,3)},...
        'MaxArea',fem.res,...
        'InGroup', group,...
        'Turns',1);
    %     mi_addblocklabel(BLKLABELSstat.xy(index,1),BLKLABELSstat.xy(index,2));
    %     mi_selectlabel(BLKLABELSstat.xy(index,1),BLKLABELSstat.xy(index,2));
    %     mi_setblockprop(BLKLABELS.materials{BLKLABELSstat.xy(index,3)}, 0, fem.res,'None', 0, group, 1);
    %     mi_clearselected;
    index=index+1;
end


end