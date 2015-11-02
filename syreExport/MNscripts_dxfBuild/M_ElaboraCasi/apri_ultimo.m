%% apri_ultimo,m - 11 06 09 - GP
%% apre ultima macchina vista e salva il percorso

%% Imposta il percorso sull'ultimo file visto
here = pwd;
if exist([here SLASH 'Ultimo.txt'],'file')
    fid=fopen([here SLASH 'Ultimo.txt'],'r');
    Percorso=fgets(fid);
    if exist(Percorso,'dir')
        here=Percorso;
    end
    fclose(fid);
    clear fid  Percorso
end

%% Apre l'ultimo motore visto
[FileName,FilePath] = uigetfile([here SLASH '*Machine.mat'],'Select the ParMachine Mat-file');
load([FilePath 'ParMachine.mat']);

%% Salva l'ultimo percorso scelto
fid=fopen([pwd SLASH 'Ultimo.txt'],'w');
for n=FilePath
    if n=='\'
        n='\\';
    end
    fprintf(fid,n);
end
fclose(fid);
clear fid n