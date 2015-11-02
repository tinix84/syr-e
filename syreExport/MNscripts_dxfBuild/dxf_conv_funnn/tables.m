function tables(fid)
%% Generazione sezione TABLES
% Scrittura intestazione TABLES
testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'SECTION';
fprintf(fid,'%s\n',testo);
testo = '2';
fprintf(fid,'%s\n',testo);
testo = 'TABLES';
fprintf(fid,'%s\n',testo);

testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'ENDSEC';
fprintf(fid,'%s\n',testo);

