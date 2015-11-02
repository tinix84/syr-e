function blocks(fid)
%% Generazione sezione BLOCKS
% Scrittura intestazione BLOCKS
testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'SECTION';
fprintf(fid,'%s\n',testo);
testo = '2';
fprintf(fid,'%s\n',testo);
testo = 'BLOCKS';
fprintf(fid,'%s\n',testo);

%scrivere testo sezione

testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'ENDSEC';
fprintf(fid,'%s\n',testo);
