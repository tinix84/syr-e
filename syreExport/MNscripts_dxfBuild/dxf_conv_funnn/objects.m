function objects(fid)
%% Generazione sezione OBJECTS
% Scrittura intestazione OBJECTS
testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'SECTION';
fprintf(fid,'%s\n',testo);
testo = '2';
fprintf(fid,'%s\n',testo);
testo = 'OBJECTS';
fprintf(fid,'%s\n',testo);

%scrivere testo sezione

testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'ENDSEC';
fprintf(fid,'%s\n',testo);

