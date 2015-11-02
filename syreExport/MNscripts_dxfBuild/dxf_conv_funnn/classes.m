function classes(fid)
%% Generazione sezione CLASSES
% Scrittura intestazione CLASSES
testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'SECTION';
fprintf(fid,'%s\n',testo);
testo = '2';
fprintf(fid,'%s\n',testo);
testo = 'CLASSES';
fprintf(fid,'%s\n',testo);

%scrivere testo sezione

testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'ENDSEC';
fprintf(fid,'%s\n',testo);

