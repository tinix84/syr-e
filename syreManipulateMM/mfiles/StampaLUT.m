% StampaLUT
fid = fopen([PATHNAME 'LUT.txt'],'w');
fprintf(fid,'Sicme SyncRel - Deflux LUT\n');
fprintf(fid,[ date '\n']);
fprintf(fid,'Nn          =  %4.0f rpm\n',w_LUT(1));
fprintf(fid,'Nmax        =  %4.0f rpm\n',w_LUT(end));

fprintf(fid,'LUT size    =  %3.0f punti\n',length(w_LUT));

fprintf(fid,'Ld_max    =  {');
for jj = 1:length(w_LUT)
    if (jj > 1)
        fprintf(fid,', %4.0f',n_Fd_max(jj));
    else
        fprintf(fid,'%4.0f',n_Fd_max(jj));
    end
end
fprintf(fid,'} \n');

fprintf(fid,'iq_max    =  {');

for jj = 1:length(w_LUT)
    if (jj > 1)
        fprintf(fid,', %4.0f',n_iq_max(jj));
    else
        fprintf(fid,'%4.0f',n_iq_max(jj));
    end
end
fprintf(fid,'} \n');
fclose(fid);