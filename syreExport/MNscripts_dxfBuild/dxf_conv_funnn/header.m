function header(fid)
%% Generazione sezione HEADER
% Scrittura intestazione HEADER
testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'SECTION';
fprintf(fid,'%s\n',testo);
testo = '  2';
fprintf(fid,'%s\n',testo);
testo = 'HEADER';
fprintf(fid,'%s\n',testo);

testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'ENDSEC';
fprintf(fid,'%s\n',testo);


