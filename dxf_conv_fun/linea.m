function linea(a,layer,color,tiplin,k,fid)
%% Scrittura dei punti che caratterizzano una linea
 
testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'LINE';
fprintf(fid,'%s\n',testo);

%definisco il Layer
testo = '8';
fprintf(fid,'%s\n',testo);
testo = layer;
fprintf(fid,'%s\n',testo);

%definisco il tipo di linea
testo = '6';
fprintf(fid,'%s\n',testo);
testo = tiplin;
fprintf(fid,'%s\n',testo);

%definisco il colore
testo = '62';
fprintf(fid,'%s\n',testo);
testo = color;
fprintf(fid,'%s\n',testo);

%definizione del primo punto
testo = '10';
fprintf(fid,'%s\n',testo);
testo = num2str(a(k,1),'%15.8f');
fprintf(fid,'%s\n',testo);
testo = '20';
fprintf(fid,'%s\n',testo);
testo = num2str(a(k,2),'%15.8f');
fprintf(fid,'%s\n',testo);
testo = '30';
fprintf(fid,'%s\n',testo);
testo = '0';
fprintf(fid,'%s\n',testo);

%definizione del secondo punto
testo = '11';
fprintf(fid,'%s\n',testo);
testo = num2str(a(k,3),'%15.8f');
fprintf(fid,'%s\n',testo);
testo = '21';
fprintf(fid,'%s\n',testo);
testo = num2str(a(k,4),'%15.8f');
fprintf(fid,'%s\n',testo);
testo = '31';
fprintf(fid,'%s\n',testo);
testo = '0';
fprintf(fid,'%s\n',testo);

