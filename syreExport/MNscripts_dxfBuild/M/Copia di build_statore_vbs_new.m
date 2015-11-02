%% VERSIONE 20 11 2011
% 07 03 08
% build_statore_vbs_new.m
% directly copied from BuildMachine1.m
% is the part that builds statore.VBS
% 27 06 08 - avvolgimento a tre strati

fid = fopen('statore.vbs','wt');

%%%%%%%%%%%%%%%%%%%%%%%%
%  INTESTAZIONE        %
%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'Call SetLocale("en-us")\n');
fprintf(fid,'%s.beginUndoGroup("Set Default Units", true)\n',CGD);
fprintf(fid,'%s.setDefaultLengthUnit("Millimeters")\n',CGD);
fprintf(fid,'%s.endUndoGroup()\n',CGD);
fprintf(fid,'%s.setCurveSmoothnessAngle (1)\n',CGD);
fprintf(fid,'\n');
fprintf(fid,'''COSTRUZIONE MEZZA CAVA\n');
fprintf(fid,'\n');

% 21 11 09 - con VDes nuovo lo statore ha un numero eccessivo di punti
rig_stat = 11;
for ii=1:rig_stat;
%     if ii == 6
%         fprintf(fid,'%s.newLine(%.15g, %.15g, %.15g, %.15g) \n',CGDGV,statore(ii,3),statore(ii,4),statore(ii,5),statore(ii,6));
%     else
        if statore(ii,col_stat)==0
            fprintf(fid,'%s.newLine(%.15g, %.15g, %.15g, %.15g) \n',CGDGV,statore(ii,1),statore(ii,2),statore(ii,3),statore(ii,4));
        else
            fprintf(fid,'%s.newArc(%.15g, %.15g, %.15g, %.15g, %.15g, %.15g) \n',CGDGV,statore(ii,1),statore(ii,2),statore(ii,3),statore(ii,4),statore(ii,5),statore(ii,6));
        end
%     end
end

fprintf(fid,'\n');
fprintf(fid,'''COSTRUZIONE STATORE\n');
fprintf(fid,'\n');

% MIRROR DI MEZZA CAVA
fprintf(fid,'%s.selectAll(infoSetSelection)\n',CGDGV);
fprintf(fid,'%s.mirrorSelectedEdges(0, 0, 1, 0, True)\n',CGDGV);

% CHIUSURA DELLA CAVA, linee laterali esterne
fprintf(fid,'%s.newLine(%.15g, %.15g, %.15g, %.15g) \n',CGDGV,statore(8,5),statore(8,6),statore(1,5),statore(1,6));
fprintf(fid,'%s.newLine(%.15g, %.15g, %.15g, %.15g) \n',CGDGV,statore(8,5),-statore(8,6),statore(1,5),-statore(1,6));

% PUNTI MEDI I MEZZA CAVA PER L'IDENTIFICAZIONE DELLE SUPERFICI DA ESTRUDERE
medS_1=(statore(9,3)+statore(10,3))/2;      % slot_air
medS_2=(statore(10,3)+statore(11,3))/2;     % slot_int
medS_3=(statore(2,1)+statore(11,3))/2;      % slot ext
medS_4=(statore(1,3)+statore(2,3))/2;       % statore: punto medio giogo
medS_5=(statore(9,3)-Mac.g*3/6);            % selezione traferro STRATO 1 (interno - tra arco a e b)
medS_6=(statore(9,3)-Mac.g/6);              % selezione traferro STRATO 2 (esterno - tra arco b e c)
%% tre strati
if (size(Magnet.avv,1) == 3)
    medS_3  =   (statore(2,1)+statore(12,3))/2;     % slot ext
    medS_7  =   (statore(11,3)+statore(12,3))/2;    % slot mid
end

% ESTRUSIONE DELLA CAVA PER L'AVVOLGIMENTO
% slot_air
fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection, Array(infoSliceSurface))\n',CGDGV,medS_1,0);
fprintf(fid,'REDIM ArrayOfValues(0)\n');
fprintf(fid,'ArrayOfValues(0)= "slot_air"\n');
fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,mat_chiu,IMCUS,IMCRV);
% slot_int
fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection, Array(infoSliceSurface))\n',CGDGV,medS_2,0);
fprintf(fid,'ArrayOfValues(0)= "slot_int"\n');
fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,mat_avv,IMCUS,IMCRV);
%% tre strati
if (size(Magnet.avv,1) == 3)
    % slot_mid
    fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection, Array(infoSliceSurface))\n',CGDGV,medS_7,0);
    fprintf(fid,'ArrayOfValues(0)= "slot_mid"\n');
    fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,mat_avv,IMCUS,IMCRV);
end
% slot_ext
fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection, Array(infoSliceSurface))\n',CGDGV,medS_3,0);
fprintf(fid,'ArrayOfValues(0)= "slot_ext"\n');
fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,mat_avv,IMCUS,IMCRV);
% statore
fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection, Array(infoSliceSurface))\n',CGDGV,medS_4,0);
fprintf(fid,'REDIM ArrayOfValues(0)\n');
fprintf(fid,'ArrayOfValues(0)= "statore"\n');
fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,mat_s,IMCUS,IMCRV);

% ROTAZIONE DELLA CAVA (occhio!! ang_asse_avv dipende da convenzione Mac.avv)
if (size(Magnet.avv,1) == 2)
    fprintf(fid,'REDIM Arrayrotate(3)\n');
elseif (size(Magnet.avv,1) == 3)
    fprintf(fid,'REDIM Arrayrotate(4)\n');
    fprintf(fid,'Arrayrotate(4)="slot_mid"\n');
end
fprintf(fid,'Arrayrotate(0)="slot_air"\n');
fprintf(fid,'Arrayrotate(1)="slot_int"\n');
fprintf(fid,'Arrayrotate(2)="slot_ext"\n');
fprintf(fid,'Arrayrotate(3)="statore"\n');
fprintf(fid,'%s.rotateComponent (Arrayrotate, 0, 0, 0, 0, 0, 1, %d, 1)\n',CGD,(-ang_asse_avv)*180/pi);


%% ruota e copia cava con nome (Copy of slot_int # -> slot_1_#)
for jj=1:1:Mac.Qs-1;
    ang_di_rotaz=ang_s*180/pi;
    fprintf(fid,'%s.copyComponent(Arrayrotate,1)\n',CGD);
    fprintf(fid,'%s.rotateComponent (Arrayrotate, 0, 0, 0, 0, 0, 1, %d, 1)\n',CGD,ang_di_rotaz);
    fprintf(fid,'%s.renameObject("Copy of slot_air #%d", "slot_air_%d")\n',CGD,jj,jj);
    fprintf(fid,'%s.renameObject("Copy of slot_int #%d", "slot_1_%d")\n',CGD,jj,jj);
    if (size(Magnet.avv,1) == 3)
        %% tre strati
        fprintf(fid,'%s.renameObject("Copy of slot_mid #%d", "slot_2_%d")\n',CGD,jj,jj);
        fprintf(fid,'%s.renameObject("Copy of slot_ext #%d", "slot_3_%d")\n',CGD,jj,jj);
    else
        %% due strati
        fprintf(fid,'%s.renameObject("Copy of slot_ext #%d", "slot_2_%d")\n',CGD,jj,jj);
    end
    fprintf(fid,'%s.renameObject("Copy of statore #%d", "statore_%d")\n',CGD,jj,jj);
    
end

%% cava master diventa l'ultima cava
fprintf(fid,'%s.renameObject("slot_air", "slot_air_%d")\n',CGD,Mac.Qs);
fprintf(fid,'%s.renameObject("slot_int", "slot_1_%d")\n',CGD,Mac.Qs);
if (size(Magnet.avv,1) == 3)
    %% tre strati
    fprintf(fid,'%s.renameObject("slot_mid", "slot_2_%d")\n',CGD,Mac.Qs);
    fprintf(fid,'%s.renameObject("slot_ext", "slot_3_%d")\n',CGD,Mac.Qs);
else
    %% due strati
    fprintf(fid,'%s.renameObject("slot_ext", "slot_2_%d")\n',CGD,Mac.Qs);
end
fprintf(fid,'%s.renameObject("statore", "statore_%d")\n',CGD,Mac.Qs);

%% traferro con doppio strato di aria
%% CALCOLO DEI PUNTI FINALI DELLO STATORE
%% allineato allo statore
ang_fine_s=ang_s*Mac.Qs;
ang_inizio_s=0;
punto1a=[(statore(9,3)-Mac.g*2/3)*cos(ang_inizio_s) (statore(9,3)-Mac.g*2/3)*sin(ang_inizio_s)];  % arco interno (a), inizio
punto2a=[(statore(9,3)-Mac.g*2/3)*cos(ang_fine_s) (statore(9,3)-Mac.g*2/3)*sin(ang_fine_s)];      % arco interno, fine
punto1b=[(statore(9,3)-Mac.g/3)*cos(ang_inizio_s) (statore(9,3)-Mac.g/3)*sin(ang_inizio_s)];      % arco intermedio (b), inizio
punto2b=[(statore(9,3)-Mac.g/3)*cos(ang_fine_s) (statore(9,3)-Mac.g/3)*sin(ang_fine_s)];          % arco intermedio, fine
punto1c=[statore(9,3)*cos(ang_inizio_s) statore(9,3)*sin(ang_inizio_s)];                          % arco esterno (c), inizio
punto2c=[statore(9,3)*cos(ang_fine_s) statore(9,3)*sin(ang_fine_s)];                              % arco esterno (c), fine

fprintf(fid,'%s.newLine(%.15g, %.15g, %.15g, %.15g) \n',CGDGV,punto1a(1),punto1a(2),punto1c(1),punto1c(2)); % linea di chiusura traferro
fprintf(fid,'%s.newLine(%.15g, %.15g, %.15g, %.15g) \n',CGDGV,punto2a(1),punto2a(2),punto2c(1),punto2c(2)); % linea di chiusura traferro
fprintf(fid,'%s.newArc(%.15g, %.15g, %.15g, %.15g, %.15g, %.15g) \n',CGDGV,statore(1,1),statore(1,2),punto1a(1),punto1a(2),punto2a(1),punto2a(2));  % arco (a)
fprintf(fid,'%s.newArc(%.15g, %.15g, %.15g, %.15g, %.15g, %.15g) \n',CGDGV,statore(1,1),statore(1,2),punto1b(1),punto1b(2),punto2b(1),punto2b(2));  % arco (b)
fprintf(fid,'%s.newArc(%.15g, %.15g, %.15g, %.15g, %.15g, %.15g) \n',CGDGV,statore(1,1),statore(1,2),punto1c(1),punto1c(2),punto2c(1),punto2c(2));  % arco (c)

% ESTRUSIONE TRAFERRO
fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection)\n',CGDGV,medS_5,0);        % selezione traferro STRATO 1 (interno - tra arco a e b)
fprintf(fid,'%s.selectAt(%.15g, %.15g, infoToggleInSelection)\n',CGDGV,medS_6,0);   % selezione traferro STRATO 2 (esterno - tra arco a e b)
fprintf(fid,'REDIM ArrayOfValues(1)\n');
fprintf(fid,'ArrayOfValues(0)= "%s"\n','air_gap_stator1');
fprintf(fid,'ArrayOfValues(1)= "%s"\n','air_gap_stator2');
fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", infoMakeComponentRemoveVertices)\n',CGDGV,Mac.l,mat_tr);

% ROTAZIONE DEL TRAFERRO DI STATORE (-ang_asse_avv-ang_s/2)
fprintf(fid,'%s.rotateComponent (ArrayOfValues, 0, 0, 0, 0, 0, 1, %d, 1)\n',CGD,-(ang_asse_avv+ang_s/2)*180/pi);

%ELIMINO LE LINEE DI COSTRUZIONE DEL TRAFERRO
fprintf(fid,'%s.selectAll(infoSetSelection, Array(infoSliceLine, infoSliceArc))\n',CGDGV);
fprintf(fid,'%s.deleteSelection()\n',CGDGV);

%MESH TRAFERRO
fprintf(fid,'REDIM ArrayOfValues(0)\n');
fprintf(fid,'ArrayOfValues(0)= "air_gap_stator1,Face#2,Edge#3"\n');
fprintf(fid,'%s.assignMeshEdgeSubdivisions(ArrayOfValues, "Type=Uniform;Subdivisions=%d;DensityRatio=1")\n',CGD,Mac.N_airgap/2);

fprintf(fid,'\n');
fprintf(fid,'''COSTRUZIONE COIL\n');
fprintf(fid,'\n');

%COSTRUZIONE COIL U
fprintf(fid,'''fase U\n');
coil_name = 'U';
coil_number = 1;
make_coil_vbs;

%COSTRUZIONE COIL V
fprintf(fid,'''fase V\n');
coil_name = 'V';
coil_number = 2;
make_coil_vbs;

%COSTRUZIONE COIL W
fprintf(fid,'''fase W\n');
coil_name = 'W';
coil_number = 3;
make_coil_vbs;

fclose(fid);