%% VERSIONE 20 11 2011
% make_coil_vbs.m   30 06 08    GMP
% make coil U V W
% called by build_statore_vbs.m


[strato_pos,cava_pos]=find(Mac.avv == coil_number);
[strato_neg,cava_neg]=find(Mac.avv == -coil_number);

%% vettore segno dei lati
u=[];
% fase U, pos
if ~isempty(strato_pos)
    for jj=1:size(strato_pos,1);
        if jj == 1
            %% selezione primo lato
            fprintf(fid,'%s.selectObject("slot_%d_%d", infoSetSelection)\n',CGDGV,strato_pos(jj),cava_pos(jj));

        else
            %% selezione altri lati
            fprintf(fid,'%s.selectObject("slot_%d_%d", %s)\n',CGDGV,strato_pos(jj),cava_pos(jj),ITIS);
        end
        u=[u ; 1];
    end
end
% fase U, neg
if ~isempty(strato_neg)
    for jj=1:size(strato_neg,1);
        if jj == 1
            %% selezione primo lato
            fprintf(fid,'%s.selectObject("slot_%d_%d", infoSetSelection)\n',CGDGV,strato_neg(jj),cava_neg(jj));
        else
            %% selezione altri lati
            fprintf(fid,'%s.selectObject("slot_%d_%d", %s)\n',CGDGV,strato_neg(jj),cava_neg(jj),ITIS);
        end
        u=[u ; -1];
    end
end

%% make coil U
%% selezione lati positivi
fprintf(fid,'REDIM ArrayOfValues(%d)\n',(size(strato_pos,1)+size(strato_neg,1)-1));
for jj=1:size(strato_pos,1);
    fprintf(fid,'ArrayOfValues(%d)= "slot_%d_%d"\n',jj-1,strato_pos(jj),cava_pos(jj));
end
%% se ci sono solo lati negativi (es. fase W)
if isempty(strato_pos)
    jj = 0;
end
%% selezione lati negativi
for kk=1:size(strato_neg,1);
    jj = jj +1; %% contatore totale
    fprintf(fid,'ArrayOfValues(%d)= "slot_%d_%d"\n',jj-1,strato_neg(kk),cava_neg(kk));
end
fprintf(fid,'%s.makeSimpleCoil(1, ArrayOfValues)\n',CGD);
fprintf(fid,'%s.renameObject("Coil#%g", "%s")\n',CGD,coil_number,coil_name);
fprintf(fid,'%s.setCoilNumberOfTurns("%s", 6)\n',CGD,coil_name);
for jj=1:size(u,1);
    if (rem(jj,2) ~= 0) && (sign(u(jj)) == -1)
        fprintf(fid,'%s.reverseCoilSide("%s", %d)\n',CGD,coil_name,jj);
    end
    if (rem(jj,2) == 0) && (sign(u(jj)) == 1)
        fprintf(fid,'%s.reverseCoilSide("%s", %d)\n',CGD,coil_name,jj);
    end
end
%% assegna proprietà coil U
fprintf(fid,'%s.selectObject("%s", infoSetSelection)\n',CGDGV,coil_name);
fprintf(fid,'%s.beginUndoGroup("Set %s Properties", true)\n',CGD,coil_name);
fprintf(fid,'REDIM ArrayOfValues(5)\n');
fprintf(fid,'ArrayOfValues(0)= 0\n');
fprintf(fid,'ArrayOfValues(1)= 33.3333\n');
fprintf(fid,'ArrayOfValues(2)= 50\n');
fprintf(fid,'ArrayOfValues(3)= 0\n');
fprintf(fid,'ArrayOfValues(4)= 0\n');
fprintf(fid,'ArrayOfValues(5)= 0\n');
fprintf(fid,'%s.setSourceWaveform("%s","SIN", ArrayOfValues)\n',CGD,coil_name);
fprintf(fid,'%s.endUndoGroup()\n',CGD);