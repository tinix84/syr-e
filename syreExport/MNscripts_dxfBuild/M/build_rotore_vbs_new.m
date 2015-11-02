%% VERSIONE 20 11 2011
% 07 03 08

% build_rotore_vbs.m
% output: rotore.VBS

% revisions:
% - 08 05 08  plasto-ferrite barriers
% - 27 06 08  shaft - opzionale
% - 21 12 09  disegna anche SMPM (aggiunta la variabile IPM)

% LEGENDA
% CODIFICA DEL TIPO DI STRATO E NUMERO DI RIGHE DELLA MATRICE rotore PER STRATO
%   siPL    = 1    Plastoferrite
%   siP     = 1    Ponticello centrale
%   siA     = 1    Aria Centrale
%   Codice per ciascuno strato: 100*SiPL + 10*SiP + 1*SiA
%   Es.                           N_comp    N Tratti (righe matrice)    N Tratti SE Ultimo
%   000     No ponticelli         1         0 + 8 + 1  = 9              0 + 6 + 1 = 7
%   001     No ponticelli, rett   1         2 + 8 + 3  = 13             2 + 6 + 3 = 11
%   010     Un ponticello         2         3 + 8 + 1  = 12             3 + 6 + 1 = 10
%   011     Due ponticelli        3         10 + 8 + 1 = 19             10 + 6 + 1 = 17
%   100,101 ..  Uguale ma con plastoferrite (due linee in più, aggiunte al fondo di ciascuno strato nella matrice rotore)

%%
% IPM o SMPM ?? (20 12 09)
IPM = isfield(Mac,'tipo_strati');
addpath E:\Matlab_functions
%%
fid = fopen('rotore.vbs','wt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INTESTAZIONE                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'Call SetLocale("en-us")\n');
fprintf(fid,'%s.beginUndoGroup("Set Default Units", true)\n',CGD);
fprintf(fid,'%s.setDefaultLengthUnit("Millimeters")\n',CGD);
fprintf(fid,'%s.endUndoGroup()\n',CGD);
fprintf(fid,'%s.setCurveSmoothnessAngle (1)\n',CGD);
fprintf(fid,'\n');
fprintf(fid,'''COSTRUZIONE ROTORE\n');
fprintf(fid,'\n');

%% DISEGNA LE BARRIERE DI ROTORE (META' SUPERIORE)
rig_rot = size(rotore,1);
for ii=1:rig_rot/2;
    if rotore(ii,col_rot)==0
        fprintf(fid,'%s.newLine(%.15g, %.15g, %.15g, %.15g) \n',CGDGV,rotore(ii,1),rotore(ii,2),rotore(ii,3),rotore(ii,4));
    elseif rotore(ii,col_rot)==1
        fprintf(fid,'%s.newArc(%.15g, %.15g, %.15g, %.15g, %.15g, %.15g) \n',CGDGV,rotore(ii,1),rotore(ii,2),rotore(ii,3),rotore(ii,4),rotore(ii,5),rotore(ii,6));
    else
        fprintf(fid,'%s.newArc(%.15g, %.15g, %.15g, %.15g, %.15g, %.15g) \n',CGDGV,rotore(ii,1),rotore(ii,2),rotore(ii,5),rotore(ii,6),rotore(ii,3),rotore(ii,4));
    end
end

if IPM
    %% DISEGNA LE BARRIERE DI ROTORE (META' INFERIORE)
    for ii=rig_rot/2+1:1:rig_rot;
        if rotore(ii,col_rot)==0
            fprintf(fid,'%s.newLine(%.15g, %.15g, %.15g, %.15g) \n',CGDGV,rotore(ii,1),rotore(ii,2),rotore(ii,3),rotore(ii,4));
        elseif rotore(ii,col_rot)==1
            fprintf(fid,'%s.newArc(%.15g, %.15g, %.15g, %.15g, %.15g, %.15g) \n',CGDGV,rotore(ii,1),rotore(ii,2),rotore(ii,3),rotore(ii,4),rotore(ii,5),rotore(ii,6));
        else
            fprintf(fid,'%s.newArc(%.15g, %.15g, %.15g, %.15g, %.15g, %.15g) \n',CGDGV,rotore(ii,1),rotore(ii,2),rotore(ii,5),rotore(ii,6),rotore(ii,3),rotore(ii,4));
        end
    end
end

if IPM
    %% TEMP !!! (19 02 10)
    fprintf(fid,'%s.selectAll(infoSetSelection, Array(infoSliceLine, infoSliceArc))\n',CGDGV);
    fprintf(fid,'%s.rotateSelectedEdges(0, 0, %d, False)\n',CGDGV,90/Mac.p);
end

%%  COMPLETA IL DISEGNO
% CALCOLA I PUNTI ESTREMI DEL ROTORE

%% mod 16 marzo 2010 - aggiunte macchine non antiperiodiche
% ang_inizio_r = 0;
if Mac.ns == Mac.Qs
    % EVEN periodic
    ang_fine_r = 2*ang_r*Mac.ps * 2;
    repetitions_r = 2*Mac.ps;
else
    % ODD periodic
    ang_fine_r = 2*ang_r*Mac.ps;
    repetitions_r = Mac.ps;
end

% arco interno ferro
if (add_shaft)
    % albero: spesso metà di raggio interno Ar
    xy_Ar_0 = 0.5*[Mac.Ar 0];
    xy_Ar_1 = 0.5*[Mac.Ar*cos(ang_fine_r) Mac.Ar*sin(ang_fine_r)];
else
    xy_Ar_0 = [rotore(1,5) rotore(1,6)];
    xy_Ar_1 = [rotore(1,5) -rotore(1,6)];
end

% arco esterno rotore
xy_xR_0=[(Mac.xR-Mac.g)*cos(0) (Mac.xR-Mac.g)*sin(0)];
xy_xR_1=[(Mac.xR-Mac.g)*cos(ang_fine_r) (Mac.xR-Mac.g)*sin(ang_fine_r)];
% arco esterno traferro (al rotore compete 1/3 di g)
xy_xR_g_0=[(Mac.xR-Mac.g*2/3)*cos(0) (Mac.xR-Mac.g*2/3)*sin(0)];
xy_xR_g_1=[(Mac.xR-Mac.g*2/3)*cos(ang_fine_r) (Mac.xR-Mac.g*2/3)*sin(ang_fine_r)];
    
    
%% DISEGNA LE linee di CHIUSURA
fprintf(fid,'%s.newLine(%.15g, %.15g, %.15g, %.15g) \n',CGDGV,xy_Ar_0(1),xy_Ar_0(2),xy_xR_0(1),xy_xR_0(2));
fprintf(fid,'%s.newLine(%.15g, %.15g, %.15g, %.15g) \n',CGDGV,xy_Ar_1(1),xy_Ar_1(2),xy_xR_1(1),xy_xR_1(2));
if (add_shaft)
    fprintf(fid,'%s.newArc(%.15g, %.15g, %.15g, %.15g, %.15g, %.15g) \n',CGDGV,0,0,xy_Ar_0(1),xy_Ar_0(2),xy_Ar_1(1),xy_Ar_1(2));
end

% air_gap_rotor
fprintf(fid,'%s.newLine(%.15g, %.15g, %.15g, %.15g) \n',CGDGV,xy_xR_0(1),xy_xR_0(2),xy_xR_g_0(1),xy_xR_g_0(2));
fprintf(fid,'%s.newLine(%.15g, %.15g, %.15g, %.15g) \n',CGDGV,xy_xR_1(1),xy_xR_1(2),xy_xR_g_1(1),xy_xR_g_1(2));
fprintf(fid,'%s.newArc(%.15g, %.15g, %.15g, %.15g, %.15g, %.15g)\n',CGDGV,rotore(1,1),rotore(1,2),xy_xR_g_0(1),xy_xR_g_0(2),xy_xR_g_1(1),xy_xR_g_1(2));

if IPM
    %% TEMP !!! (19 02 10)
    fprintf(fid,'%s.selectAll(infoSetSelection, Array(infoSliceLine, infoSliceArc))\n',CGDGV);
    fprintf(fid,'%s.rotateSelectedEdges(0, 0, %d, False)\n',CGDGV,-90/Mac.p);
end

%%%
%% ESTRUSIONI
%%%
fprintf(fid,'\n');
fprintf(fid,'''ESTRUSIONE ROTORE\n');
fprintf(fid,'\n');

medR_r = Mac.Ar + 0.1;
angle_pole = ang_r + (0:repetitions_r-1) * (2 * ang_r);
[medR_x,medR_y] = pol2cart(angle_pole,medR_r);

% ESTRUSIONE delle BARRIERE/MAGNETI
if IPM
    a=1;m=1;    % indice lati aria, indice lati magneti
    for jj=1:1:n_strati
        switch Mac.tipo_strati(jj)
            case {0,1}   % 000 strato no ponticello (1 segmento)
                if jj == n_strati    % 7 lati
                    med=[(rotore(1+a,1)+rotore(7+a,3))/2 (rotore(1+a,2)+rotore(7+a,4))/2];
                    fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection, Array(infoSliceSurface))\n',CGDGV,med(1),med(2));
                    fprintf(fid,'REDIM ArrayOfValues(0)\n');
                    fprintf(fid,'ArrayOfValues(0)= "barrier_%d"\n',jj);
                    fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,mat_barr,IMCUS,IMCRV);
                else                 % 9 lati
                    med=[(rotore(1+a,1)+rotore(9+a,3))/2 (rotore(1+a,2)+rotore(9+a,4))/2];
                    fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection, Array(infoSliceSurface))\n',CGDGV,med(1),med(2));
                    fprintf(fid,'REDIM ArrayOfValues(0)\n');
                    fprintf(fid,'ArrayOfValues(0)= "barrier_%d"\n',jj);
                    fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,mat_barr,IMCUS,IMCRV);
                end
                a=a+9;
            case 10   % 010 strato con un ponticello (2 segmenti)
                fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection, Array(infoSliceSurface))\n',CGDGV,rotore(1+a,1),rotore(1+a,2));
                fprintf(fid,'REDIM ArrayOfValues(0)\n');
                fprintf(fid,'ArrayOfValues(0)= "barrier_%d_a"\n',jj);
                fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,mat_barr,IMCUS,IMCRV);
                fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection, Array(infoSliceSurface))\n',CGDGV,rotore(1+a,1),-rotore(1+a,2));
                fprintf(fid,'REDIM ArrayOfValues(0)\n');
                fprintf(fid,'ArrayOfValues(0)= "barrier_%d_b"\n',jj);
                fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,mat_barr,IMCUS,IMCRV);
                if jj == n_strati    % 10 lati
                    a=a+10;
                else                 % 12 lati
                    a=a+12;
                end
            case 11   % 011 strato con due ponticelli (3 segmenti)
                
                med_2=[(rotore(1+a,1)+rotore(7+a,3))/2 (rotore(1+a,2)+rotore(7+a,4))/2];
                fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection, Array(infoSliceSurface))\n',CGDGV,rotore(8+a,1),rotore(8+a,2));
                fprintf(fid,'REDIM ArrayOfValues(0)\n');
                fprintf(fid,'ArrayOfValues(0)= "barrier_%d_a"\n',jj);
                fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,mat_barr,IMCUS,IMCRV);
                fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection, Array(infoSliceSurface))\n',CGDGV,med_2(1),med_2(2));
                fprintf(fid,'REDIM ArrayOfValues(0)\n');
                fprintf(fid,'ArrayOfValues(0)= "barrier_%d_b"\n',jj);
                fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,mat_barr,IMCUS,IMCRV);
                fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection, Array(infoSliceSurface))\n',CGDGV,rotore(8+a,1),-rotore(8+a,2));
                fprintf(fid,'REDIM ArrayOfValues(0)\n');
                fprintf(fid,'ArrayOfValues(0)= "barrier_%d_c"\n',jj);
                fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,mat_barr,IMCUS,IMCRV);
                
                if jj == n_strati        % 17 lati
                    a=a+17;
                else                     % 19 lati
                    a=a+19;
                end
            case 100   % 100 plastoferrite, no ponticello (1 segmento)
                if jj == n_strati    % 7 lati - ultimo strato
                    med=[(rotore(1+a,1)+rotore(7+a,3))/2 (rotore(1+a,2)+rotore(7+a,4))/2];
                    fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection,Array(infoSliceSurface))\n',CGDGV,med(1),med(2));  % punto medio lato di base
                    fprintf(fid,'%s.selectAt(%.15g, %.15g, infoToggleInSelection,Array(infoSliceSurface))\n',CGDGV,mean(rotore(6+a,3:2:5)), ...
                        mean(rotore(6+a,4:2:6)));                                                                             % media coordinate arco esterno (lato 6)
                    fprintf(fid,'%s.selectAt(%.15g, %.15g, infoToggleInSelection,Array(infoSliceSurface))\n',CGDGV, ...
                        mean([rotore(3+a,1:2:3),rotore(8+a,1:2:3),rotore(2+a,3:2:5),rotore(4+a,3:2:5)]), ...
                        mean([rotore(3+a,2:2:4),rotore(8+a,2:2:4),rotore(2+a,4:2:6),rotore(4+a,4:2:6)]));                       % baricentro del trapezio punta banana (lati 3,8, e archi 2,4)
                    fprintf(fid,'%s.selectAt(%.15g, %.15g, infoToggleInSelection,Array(infoSliceSurface))\n',CGDGV,mean(rotore(6+a,3:2:5)), ...
                        -mean(rotore(6+a,4:2:6)));                                                                            % media coordinate arco esterno (lato 6 negato)
                    fprintf(fid,'%s.selectAt(%.15g, %.15g, infoToggleInSelection,Array(infoSliceSurface))\n',CGDGV, ...
                        mean([rotore(3+a,1:2:3),rotore(8+a,1:2:3),rotore(2+a,3:2:5),rotore(4+a,3:2:5)]), ...
                        -mean([rotore(3+a,2:2:4),rotore(8+a,2:2:4),rotore(2+a,4:2:6),rotore(4+a,4:2:6)]));                      % baricentro del trapezio punta banana (lati 3,8, e archi 2,4)                    fprintf(fid,'REDIM ArrayOfValues(4)\n');
                    fprintf(fid,'REDIM ArrayOfValues(4)\n');
                    fprintf(fid,'ArrayOfValues(0)= "plasto_%d0"\n',jj);
                    fprintf(fid,'ArrayOfValues(1)= "plasto_%d1_a"\n',jj);
                    fprintf(fid,'ArrayOfValues(2)= "plasto_%d2_a"\n',jj);
                    fprintf(fid,'ArrayOfValues(3)= "plasto_%d1_b"\n',jj);
                    fprintf(fid,'ArrayOfValues(4)= "plasto_%d2_b"\n',jj);
                    fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s;Type=Uniform;Direction=[%g,%g,0]", infoMakeComponentRemoveVertices)\n',CGDGV,Mac.l,mat_barr,1,0);
                    for ii = 1:(Mac.N_comp_strati(jj))
                        fprintf(fid,'%s.assignMaterial(ArrayOfValues(%g), "Name=%s;Type=Uniform;Direction=[%g,%g,0]")\n',CGD,ii-1,mat_barr,real(PL_MxMy(jj,ii)),imag(PL_MxMy(jj,ii)));
                    end
                else                 %9 lati
                    med=[(rotore(1+a,1)+rotore(9+a,3))/2 (rotore(1+a,2)+rotore(9+a,4))/2];
                    fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection,Array(infoSliceSurface))\n',CGDGV,med(1),med(2));  % punto medio lato di base
                    med1(1)=mean([rotore(9+a,1) rotore(1+a,3) rotore(7+a,3)]);
                    med1(2)=mean([rotore(9+a,2) rotore(1+a,4) rotore(7+a,4)]);
                    fprintf(fid,'%s.selectAt(%.15g, %.15g, infoToggleInSelection,Array(infoSliceSurface))\n',CGDGV,med1(1),med1(2));
                    med2(1)=mean([rotore(3+a,1) rotore(7+a,1) rotore(7+a,3)]);
                    med2(2)=mean([rotore(3+a,2) rotore(7+a,2) rotore(7+a,4)]);
                    fprintf(fid,'%s.selectAt(%.15g, %.15g, infoToggleInSelection,Array(infoSliceSurface))\n',CGDGV,med2(1),med2(2));
                    fprintf(fid,'%s.selectAt(%.15g, %.15g, infoToggleInSelection,Array(infoSliceSurface))\n',CGDGV,med1(1),-med1(2));
                    fprintf(fid,'%s.selectAt(%.15g, %.15g, infoToggleInSelection,Array(infoSliceSurface))\n',CGDGV,med2(1),-med2(2));
                    fprintf(fid,'REDIM ArrayOfValues(4)\n');
                    fprintf(fid,'ArrayOfValues(0)= "plasto_%d0"\n',jj);
                    fprintf(fid,'ArrayOfValues(1)= "plasto_%d1_a"\n',jj);
                    fprintf(fid,'ArrayOfValues(2)= "plasto_%d2_a"\n',jj);
                    fprintf(fid,'ArrayOfValues(3)= "plasto_%d1_b"\n',jj);
                    fprintf(fid,'ArrayOfValues(4)= "plasto_%d2_b"\n',jj);
                    fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s;Type=Uniform;Direction=[%g,%g,0]", infoMakeComponentRemoveVertices)\n',CGDGV,Mac.l,mat_barr,1,0);
                    for ii = 1:(Mac.N_comp_strati(jj))
                        fprintf(fid,'%s.assignMaterial(ArrayOfValues(%g), "Name=%s;Type=Uniform;Direction=[%g,%g,0]")\n',CGD,ii-1,mat_barr,real(PL_MxMy(jj,ii)),imag(PL_MxMy(jj,ii)));
                    end
                    a=a+9+2;    % +2 x plastoferrite
                end
            case 110   % 110 ponticello (2 segmenti), plastoferrite
                fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection,Array(infoSliceSurface))\n',CGDGV,rotore(1+a,1),rotore(1+a,2));      % centro raccordo ponticello
                fprintf(fid,'%s.selectAt(%.15g, %.15g, infoToggleInSelection,Array(infoSliceSurface))\n',CGDGV,mean(rotore(11+a,3:2:5)), ...
                    mean(rotore(11+a,4:2:6)));                                                                     % media coordinate arco esterno (lato 12)
                fprintf(fid,'%s.selectAt(%.15g, %.15g, infoToggleInSelection,Array(infoSliceSurface))\n',CGDGV,rotore(9+a,1),rotore(9+a,2)); % centro raccordo superiore (lato 10)
                fprintf(fid,'%s.selectAt(%.15g, %.15g, infoToggleInSelection,Array(infoSliceSurface))\n',CGDGV,rotore(1+a,1),-rotore(1+a,2));
                fprintf(fid,'%s.selectAt(%.15g, %.15g, infoToggleInSelection,Array(infoSliceSurface))\n',CGDGV,mean(rotore(11+a,3:2:5)),-mean(rotore(11+a,4:2:6)));   % media coordinate lato 12 (arco)
                fprintf(fid,'%s.selectAt(%.15g, %.15g, infoToggleInSelection,Array(infoSliceSurface))\n',CGDGV,rotore(9+a,1),-rotore(9+a,2)); % centro lato 10 (arco)
                fprintf(fid,'REDIM ArrayOfValues(5)\n');
                fprintf(fid,'ArrayOfValues(0)= "plasto_%d0_a"\n',jj);
                fprintf(fid,'ArrayOfValues(1)= "plasto_%d1_a"\n',jj);
                fprintf(fid,'ArrayOfValues(2)= "plasto_%d2_a"\n',jj);
                fprintf(fid,'ArrayOfValues(3)= "plasto_%d0_b"\n',jj);
                fprintf(fid,'ArrayOfValues(4)= "plasto_%d1_b"\n',jj);
                fprintf(fid,'ArrayOfValues(5)= "plasto_%d2_b"\n',jj);
                fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s;Type=Uniform;Direction=[%g,%g,0]", infoMakeComponentRemoveVertices)\n',CGDGV,Mac.l,mat_barr,1,0);
                for ii = 1:(Mac.N_comp_strati(jj))
                    fprintf(fid,'%s.assignMaterial(ArrayOfValues(%g), "Name=%s;Type=Uniform;Direction=[%g,%g,0]")\n',CGD,ii-1,mat_barr,real(PL_MxMy(jj,ii)),imag(PL_MxMy(jj,ii)));
                end
                
                if jj == n_strati    % ultimo strato: 10 lati +2 x plastoferrite
                    a=a+10+2;
                else                 % 12 lati +2 x plastoferrite
                    a=a+12+2;
                end
            otherwise
                msgbox('Error!!! blurp0');
        end
    end
end

%% DISEGNA I MAGNETI
fprintf(fid,'\n');
fprintf(fid,'''DISEGNA I MAGNETI\n');
fprintf(fid,'\n');

if ~isempty(Mac.magneti)
    magneti = Mac.magneti;    
    for ii=1:size(magneti,1)
        if magneti(ii,col_rot)==0
            fprintf(fid,'%s.newLine(%.15g, %.15g, %.15g, %.15g) \n',CGDGV,magneti(ii,1),magneti(ii,2),magneti(ii,3),magneti(ii,4));
        else
            fprintf(fid,'%s.newArc(%.15g, %.15g, %.15g, %.15g, %.15g, %.15g) \n',CGDGV,magneti(ii,1),magneti(ii,2),magneti(ii,3),magneti(ii,4),magneti(ii,5),magneti(ii,6));
        end
%         if IPM
%             fprintf(fid,'%s.newLine(%.15g, %.15g, %.15g, %.15g) \n',CGDGV,magneti(ii,1),-magneti(ii,2),magneti(ii,3),-magneti(ii,4));
%         end
    end
end

if ~IPM
    % aria alleggerimento rotore
    linee_aria_rot = rotore([5:8],:);
    punti_aria_rot = [linee_aria_rot(1:2,3:6);linee_aria_rot(3:4,1:4)];
    
    aria_rot_punto_ext1 = linee_aria_rot(1,3:4);
    aria_rot_punto_ext2 = linee_aria_rot(1,5:6);
    [temp_a_ext1,temp_r_ext] = cart2pol(aria_rot_punto_ext1(1),aria_rot_punto_ext1(2));
    [temp_a_ext2,temp_r_ext] = cart2pol(aria_rot_punto_ext2(1),aria_rot_punto_ext2(2));
    temp_a = mean([temp_a_ext1 temp_a_ext2]);
    aria_rot_punto_int1 = linee_aria_rot(2,3:4);
    aria_rot_punto_int2 = linee_aria_rot(2,5:6);
    [temp,temp_r_int] = cart2pol(aria_rot_punto_int1(1),aria_rot_punto_int1(2));
    [x_aria_rot,y_aria_rot] = pol2cart(temp_a,mean([temp_r_int temp_r_ext]))
    
    [temp_x,temp_y] = rot_point(x_aria_rot,y_aria_rot,angle_pole-angle_pole(1));
    fprintf(fid,'REDIM ArrayOfValues(0)\n');
    for jj = 1:length(angle_pole)
        fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection)\n',CGDGV,temp_x(jj),temp_y(jj));
        sel_name = ['air_rot_' num2str(jj)];
        fprintf(fid,'ArrayOfValues(0)= "%s"\n',sel_name);
        fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,'AIR',IMCUS,IMCRV);
    end
    % ARIA INTERPOLO (19 02 10)
    temp_r = mean([rotore(9,1) rotore(9,3)]);
    [temp_x1,temp_y1] = pol2cart(angle_pole-angle_pole(1)*0.95,temp_r);
    [temp_x2,temp_y2] = pol2cart(angle_pole+angle_pole(1)*0.95,temp_r);
    for jj = 1:length(angle_pole)
        % aria a
        fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection, Array(infoSliceSurface))\n',CGDGV,temp_x1(jj),temp_y1(jj));
        fprintf(fid,'ArrayOfValues(0)= "%s"\n',['air_rot_' num2str(jj) '_a']);
        fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,'AIR',IMCUS,IMCRV);
        % aria b
        fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection, Array(infoSliceSurface))\n',CGDGV,temp_x2(jj),temp_y2(jj));
        fprintf(fid,'ArrayOfValues(0)= "%s"\n',['air_rot_' num2str(jj) '_b']);
        fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,'AIR',IMCUS,IMCRV);
    end
end

% CALCOLA I PUNTI DI ESTRUSIONE DEL ROTORE (parte in ferro)
medR_tr = Mac.xR - 5/6*Mac.g;
if IPM
    medR_r  =(rotore(1,3)+ 1);
else
    medR_r  = Mac.Ar + 1;
end
medR_shaft = Mac.Ar * 0.75;

% ESTRUDE IL FERRO ROTORE
fprintf(fid,'REDIM ArrayOfValues(%g)\n',length(angle_pole)-1);
for jj = 1:length(angle_pole)
    if jj == 1
        sel_option = 'infoSetSelection';
        sel_name = 'rotor';
    else
        sel_option = 'infoToggleInSelection';
        sel_name = ['Component#' num2str(jj)];
    end
    fprintf(fid,'%s.selectAt(%.15g, %.15g, %s)\n',CGDGV,medR_x(jj),medR_y(jj),sel_option);
    fprintf(fid,'ArrayOfValues(%g)= "%s"\n',jj-1,sel_name);
end
fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,mat_r,IMCUS,IMCRV);

if (add_shaft)
    % ESTRUDE ALBERO
    [x_temp,y_temp] = rot_point(medR_shaft,0,0.9*ang_r);
    fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection, Array(infoSliceSurface))\n',CGDGV,x_temp,y_temp);
    fprintf(fid,'REDIM ArrayOfValues(0)\n');
    fprintf(fid,'ArrayOfValues(0)= "%s"\n','shaft');
    fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,mat_shaft,IMCUS,IMCRV);
end

%% ESTRUDI I MAGNETI
fprintf(fid,'\n');
fprintf(fid,'''ESTRUSIONE DEI MAGNETI\n');
fprintf(fid,'\n');
if ~isempty(Mac.magneti)
    if IPM
        n_el_m = 6;
    else
        n_el_m = 4;
    end
    n_magneti_polo = size(magneti,1)/n_el_m/repetitions_r;
    % m = 1 : numero di magneti (singolo tegolino)
    m=1;
    % mp = 1 : n_magneti_polo (cambia il verso di magnetizzazione quando
    % cambia il polo)
    mp = 1;
    verso_mag = -1;
    
    for ii=1:size(magneti,1)/n_el_m;
        if IPM
            % temporaneo, vale solo qui .. per uniformare la procedure
            % magneti IPM deve essere spedito da VDes come per SMPM (4
            % punti e pedalare e non 6 ordinati a cazzo)
            n_el_m = 3; 
            baricentro_mag=[(magneti(m,1)+magneti(2+m,1))/2 0];
        else
            xy_vertici = [magneti(m,1:2:3) magneti(m+2,1:2:3)] +1i * [magneti(m,2:2:4) magneti(m+2,2:2:4)];
            r_vertici = abs(xy_vertici);
            t_vertici = angle(xy_vertici);
            r_baricentro = mean([mean(r_vertici(1:2)) mean(r_vertici(3:4))]);
            t_baricentro = mean([mean(t_vertici(1:2)) mean(t_vertici(3:4))]);
            temp = r_baricentro*exp(1i*t_baricentro);
            baricentro_mag = [real(temp) imag(temp)];
            mag_x = cos(t_baricentro); mag_y = sin(t_baricentro);
        end
        fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection, Array(infoSliceSurface))\n',CGDGV,baricentro_mag(1,1),baricentro_mag(1,2));
        fprintf(fid,'REDIM ArrayOfValues(0)\n');
        fprintf(fid,'ArrayOfValues(0)= "magnet_%d"\n',ii);
        
        if IPM
            string_magnet = [mat_magn ';Type=Uniform;Direction=[-1,0,0]'];
        else
            string_magnet = [mat_magn ';Type=Uniform;Direction=[' num2str(verso_mag * mag_x) ',' num2str(verso_mag * mag_y) ',0]'];
        end
        
        fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,string_magnet,IMCUS,IMCRV);
        %         fprintf(fid,'%s.assignMaterial(ArrayOfValues, "Name=%s", %s Or %s)\n',CGD,mat_magn,IMCUS,IMCRV);
        m = m + n_el_m;
        mp = mp +1;
        if mp > n_magneti_polo
            mp = 1;
            verso_mag = - verso_mag;
        end
    end
    if IPM
        % temporaneo, da rimuovere
        n_el_m = 6;
    end
end

% ESTRUDE IL TRAFERRO ROTORE
fprintf(fid,'%s.selectAt(%.15g, %.15g, infoSetSelection, Array(infoSliceSurface))\n',CGDGV,medR_tr,0);
fprintf(fid,'REDIM ArrayOfValues(0)\n');
fprintf(fid,'ArrayOfValues(0)= "%s"\n','air_gap_rotor');
fprintf(fid,'%s.makeComponentInALine(%g, ArrayOfValues, "Name=%s", %s Or %s)\n',CGDGV,Mac.l,mat_tr,IMCUS,IMCRV);

%     % ELIMINO LE LINEE DI COSTRUZIONE DEL ROTORE
%     fprintf(fid,'%s.selectAll(infoSetSelection, Array(infoSliceLine, infoSliceArc))\n',CGDGV);
%     fprintf(fid,'%s.deleteSelection()\n',CGDGV);

%% ROTAZIONE DEL ROTORE E COSTRUZIONE MOTION
fprintf(fid,'\n');
fprintf(fid,'''COSTRUZIONE MOTION\n');
fprintf(fid,'\n');

if IPM
    % IPM: Gli elementi da estrudere sono:
    % rotor, shaft (eventuale), 1 air_gap, sum(Mac.N_comp_strati) numero di
    % pezzi di aria o plastoferrite, size(magneti,1)/3 numero di magneti
    % -1        per il conteggio in C
    if ~isempty(Mac.magneti)
        % con magneti (IPM)
        N_comp_estr=1+add_shaft+1+sum(Mac.N_comp_strati)+size(magneti,1)/n_el_m-1;
    else
        % senza magneti (SyR)
        N_comp_estr=1+add_shaft+1+sum(Mac.N_comp_strati)-1;
    end
else
    % Gli elementi da estrudere sono:
    % rotor, shaft (eventuale), air_gap, size(magneti,1)/4 * magneti,
    % Mac.ps * air_rot, 2 * Mac.ps * air_rot_ab
    % -1        per il conteggio in C
    N_comp_estr=1+add_shaft+1+size(magneti,1)/n_el_m+repetitions_r+2*repetitions_r-1;
end
fprintf(fid,'REDIM ArrayMotion(%d)\n',N_comp_estr);
ind_arr=1;
fprintf(fid,'%s.selectObject("rotor", infoSetSelection)\n',CGDGV);
fprintf(fid,'ArrayMotion(0)= "rotor"\n');
if IPM
    for jj=1:1:n_strati;
        switch Mac.tipo_strati(jj)
            case {0,1}
                fprintf(fid,'%s.selectObject("barrier_%d",%s)\n',CGDGV,jj,ITIS);
                fprintf(fid,'ArrayMotion(%d)= "barrier_%d"\n',ind_arr,jj);
                ind_arr=ind_arr+1;
                if Mac.tipo_strati(jj) == 1 % aggiungo magnete alla motion
                    fprintf(fid,'%s.selectObject("magnet_%d",%s)\n',CGDGV,jj,ITIS);
                    fprintf(fid,'ArrayMotion(%d)= "magnet_%d"\n',ind_arr,jj);
                    ind_arr=ind_arr+1;
                end
            case 10
                fprintf(fid,'%s.selectObject("barrier_%d_a",%s)\n',CGDGV,jj,ITIS);
                fprintf(fid,'%s.selectObject("barrier_%d_b",%s)\n',CGDGV,jj,ITIS);
                
                fprintf(fid,'ArrayMotion(%d)= "barrier_%d_a"\n',ind_arr,jj);
                fprintf(fid,'ArrayMotion(%d)= "barrier_%d_b"\n',ind_arr+1,jj);
                ind_arr=ind_arr+2;
            case 11
                fprintf(fid,'%s.selectObject("barrier_%d_a",%s)\n',CGDGV,jj,ITIS);
                fprintf(fid,'%s.selectObject("barrier_%d_b",%s)\n',CGDGV,jj,ITIS);
                fprintf(fid,'%s.selectObject("barrier_%d_c",%s)\n',CGDGV,jj,ITIS);
                
                fprintf(fid,'ArrayMotion(%d)= "barrier_%d_a"\n',ind_arr,jj);
                fprintf(fid,'ArrayMotion(%d)= "barrier_%d_b"\n',ind_arr+1,jj);
                fprintf(fid,'ArrayMotion(%d)= "barrier_%d_c"\n',ind_arr+2,jj);
                if isempty(Mac.magneti) == 0    % Mac.magneti (IPM)
                    fprintf(fid,'%s.selectObject("magnet_%d",%s)\n',CGDGV,jj,ITIS);
                    fprintf(fid,'ArrayMotion(%d)= "magnet_%d"\n',ind_arr+3,jj);
                    ind_arr=ind_arr+1;
                end
                ind_arr=ind_arr+3;
            case {100,110}
                if (Mac.tipo_strati(jj) == 100)
                    fprintf(fid,'%s.selectObject("plasto_%d0",%s)\n',CGDGV,jj,ITIS);
                    fprintf(fid,'ArrayMotion(%d)= "plasto_%d0"\n',ind_arr,jj);
                else
                    fprintf(fid,'%s.selectObject("plasto_%d0_a",%s)\n',CGDGV,jj,ITIS);
                    fprintf(fid,'%s.selectObject("plasto_%d0_b",%s)\n',CGDGV,jj,ITIS);
                    fprintf(fid,'ArrayMotion(%d)= "plasto_%d0_a"\n',ind_arr,jj);
                    fprintf(fid,'ArrayMotion(%d)= "plasto_%d0_b"\n',ind_arr+1,jj);
                    ind_arr = ind_arr+1;
                    
                end
                fprintf(fid,'%s.selectObject("plasto_%d1_a",%s)\n',CGDGV,jj,ITIS);
                fprintf(fid,'%s.selectObject("plasto_%d2_a",%s)\n',CGDGV,jj,ITIS);
                fprintf(fid,'%s.selectObject("plasto_%d1_b",%s)\n',CGDGV,jj,ITIS);
                fprintf(fid,'%s.selectObject("plasto_%d2_b",%s)\n',CGDGV,jj,ITIS);
                
                fprintf(fid,'ArrayMotion(%d)= "plasto_%d1_a"\n',ind_arr+1,jj);
                fprintf(fid,'ArrayMotion(%d)= "plasto_%d2_a"\n',ind_arr+2,jj);
                fprintf(fid,'ArrayMotion(%d)= "plasto_%d1_b"\n',ind_arr+3,jj);
                fprintf(fid,'ArrayMotion(%d)= "plasto_%d2_b"\n',ind_arr+4,jj);
                ind_arr=ind_arr+5;
            otherwise
                msgbox('Error!!! blurp1');
        end         % end switch
    end                 % end for
else
    for jj = 1:size(magneti,1)/4
        fprintf(fid,'%s.selectObject("magnet_%d",%s)\n',CGDGV,jj,ITIS);
        fprintf(fid,'ArrayMotion(%d)= "magnet_%d"\n',ind_arr,jj);
        ind_arr = ind_arr + 1;
    end
    for jj = 1:repetitions_r
        fprintf(fid,'%s.selectObject("air_rot_%d",%s)\n',CGDGV,jj,ITIS);
        fprintf(fid,'ArrayMotion(%d)= "air_rot_%d"\n',ind_arr,jj);
        fprintf(fid,'%s.selectObject("air_rot_%d_a",%s)\n',CGDGV,jj,ITIS);
        fprintf(fid,'ArrayMotion(%d)= "air_rot_%d_a"\n',ind_arr+1,jj);
        fprintf(fid,'%s.selectObject("air_rot_%d_b",%s)\n',CGDGV,jj,ITIS);
        fprintf(fid,'ArrayMotion(%d)= "air_rot_%d_b"\n',ind_arr+2,jj);
        ind_arr = ind_arr + 3;
    end
end
fprintf(fid,'%s.selectObject("air_gap_rotor", %s)\n',CGDGV,ITIS);
fprintf(fid,'ArrayMotion(%d)= "air_gap_rotor"\n',ind_arr);

if (add_shaft)
    fprintf(fid,'%s.selectObject("shaft", %s)\n',CGDGV,ITIS);
    fprintf(fid,'ArrayMotion(%d)= "shaft"\n',ind_arr+1);
end

if IPM
    fprintf(fid,'%s.rotateComponent (ArrayMotion, 0, 0, 0, 0, 0, 1, %d, 1)\n',CGD,ang_r*180/pi);
end
fprintf(fid,'%s.makeMotionComponent(ArrayMotion)\n',CGD);
fprintf(fid,'%s.unselectAll()\n',CGDGV);

%PARAMETRI DELLA MOTION
fprintf(fid,'%s.selectObject("Motion#1", infoSetSelection)\n',CGDGV);
fprintf(fid,'%s.setMotionSourceType("Motion#1", infoVelocityDriven)\n',CGD);
fprintf(fid,'%s.beginUndoGroup("Set Motion#1 Properties", true)\n',CGD);
fprintf(fid,'%s.setMotionPositionAtStartup("Motion#1", 0)\n',CGD);
fprintf(fid,'%s.setMotionSpeedAtStartup("Motion#1", 0)\n',CGD);
fprintf(fid,'REDIM ArrayOfValues1(1)\n');
fprintf(fid,'ArrayOfValues1(0)= 0\n');
fprintf(fid,'ArrayOfValues1(1)= 10\n');
fprintf(fid,'REDIM ArrayOfValues2(1)\n');
fprintf(fid,'ArrayOfValues2(0)= 0\n');
fprintf(fid,'ArrayOfValues2(1)= 0\n');
fprintf(fid,'%s.setMotionPositionVsTime("Motion#1", ArrayOfValues1, ArrayOfValues2)\n',CGD);
fprintf(fid,'%s.setMotionRotaryCenter("Motion#1", Array(0, 0, 0))\n',CGD);
fprintf(fid,'%s.setMotionRotaryAxis("Motion#1", Array(0, 0, 1))\n',CGD);
fprintf(fid,'%s.endUndoGroup()\n',CGD);

%% CONDIZIONI AL CONTORNO
fprintf(fid,'\n');
fprintf(fid,'''CONDIZIONI AL CONTORNO\n');
fprintf(fid,'\n');
% 
% Boundary Condition #1 - Odd Periodicity (or Even)

% default = IPM
face_stator    = 3;
face_airgaps2  = 4;
face_airgaps1  = 4;
face_airgaps11 = 3;
face_airgapr   = 3;
face_rotor     = 3;
face_shaft     = 3;

if ~IPM
    face_airrot    = 4;
%     if round(Mac.q) == Mac.q
    if ang_fine_r <= pi/2
        % cave intere
        face_airgapr   = 4;
        face_rotor     = 4;
        face_shaft     = 4;
    else
        % cave frazionarie
        face_airgaps2  = 5;
        face_airgaps1  = 5;
        face_airgaps11 = 4;
        face_airgapr   = 5;
        face_rotor     = 5;
        face_shaft     = 5;
    end
end

if IPM
    n_elements = 5;
else
    n_elements = 6;
end
if (add_shaft)
    n_elements = n_elements + 1;
end

fprintf(fid,'REDIM ArrayOfValues(%d)\n',n_elements);

if (add_shaft)
    fprintf(fid,'%s.selectObject("shaft,Face#%d", infoSetSelection)\n',CGDGV,face_shaft);
    fprintf(fid,'ArrayOfValues(%d)= "shaft,Face#%d"\n',n_elements,face_shaft);
end
fprintf(fid,'%s.selectObject("rotor,Face#%d", infoSetSelection)\n',CGDGV,face_rotor);
fprintf(fid,'ArrayOfValues(%d)= "rotor,Face#%d"\n',n_elements-1,face_rotor);
fprintf(fid,'%s.selectObject("air_gap_rotor,Face#%d", %s)\n',CGDGV,face_airgapr,ITIS);
fprintf(fid,'ArrayOfValues(%d)= "air_gap_rotor,Face#%d"\n',n_elements-2,face_airgapr);
fprintf(fid,'%s.selectObject("air_gap_stator1,Face#%d", %s)\n',CGDGV,face_airgaps11,ITIS);
fprintf(fid,'ArrayOfValues(%d)= "air_gap_stator1,Face#%d"\n',n_elements-3,face_airgaps11);
fprintf(fid,'%s.selectObject("air_gap_stator1,Face#%d", %s)\n',CGDGV,face_airgaps1,ITIS);
fprintf(fid,'ArrayOfValues(%d)= "air_gap_stator1,Face#%d"\n',n_elements-4,face_airgaps1);
fprintf(fid,'%s.selectObject("air_gap_stator2,Face#%d", %s)\n',CGDGV,face_airgaps2,ITIS);
fprintf(fid,'ArrayOfValues(%d)= "air_gap_stator2,Face#%d"\n',n_elements-5,face_airgaps2);
fprintf(fid,'%s.selectObject("statore_1,Face#%d", %s)\n',CGDGV,face_stator,ITIS);
fprintf(fid,'ArrayOfValues(%d)= "statore_1,Face#3"\n',n_elements-6);

if ~IPM
    fprintf(fid,'%s.selectObject("air_rot_1_a,Face#%d", %s)\n',CGDGV,face_airrot,ITIS);
    fprintf(fid,'ArrayOfValues(0)= "air_rot_1_a,Face#%d"\n',face_airrot);
end

fprintf(fid,'%s.beginUndoGroup("Assign Boundary Condition")\n',CGD);
fprintf(fid,'%s.createBoundaryCondition(ArrayOfValues, "BoundaryCondition#1")\n',CGD);
fprintf(fid,'REDIM RotationAxis(2)\n');
fprintf(fid,'RotationAxis(0)= 0\n');
fprintf(fid,'RotationAxis(1)= 0\n');
fprintf(fid,'RotationAxis(2)= 1\n');
fprintf(fid,'REDIM Center(2)\n');
fprintf(fid,'Center(0)= 0\n');
fprintf(fid,'Center(1)= 0\n');
fprintf(fid,'Center(2)= 0\n');
if Mac.ns == Mac.Qs
    % EVEN periodic
    fprintf(fid,'%s.setEvenPeriodic("BoundaryCondition#1", 1, %d, RotationAxis, Null, Null, Center)\n',CGD,ang_fine_r * 180/pi);
else
    % ODD periodic
    fprintf(fid,'%s.setOddPeriodic("BoundaryCondition#1", 1, %d, RotationAxis, Null, Null, Center)\n',CGD,ang_fine_r * 180/pi);
end
fprintf(fid,'%s.endUndoGroup()\n',CGD);

% Boundary Condition #2 - Flux Tangential - Stator
% numero cave di statore
% if round(Mac.Qs) == Mac.Qs
%     repetitions = Mac.Qs - 1;
% else
%     repetitions = Mac.Qs * 2 -1;
% end
repetitions_s = Mac.Qs - 1;
for jj=1:repetitions_s
    if jj == 1
        fprintf(fid,'%s.selectObject("statore_1,Face#4", infoSetSelection)\n',CGDGV);
    else
        fprintf(fid,'%s.selectObject("statore_%d,Face#4", %s)\n',CGDGV,jj,ITIS);
    end
end
fprintf(fid,'%s.beginUndoGroup("Assign Boundary Condition")\n',CGD);
fprintf(fid,'REDIM ArrayOfValues(%d)\n',Mac.Qs-1);
for jj=1:Mac.Qs
    fprintf(fid,'ArrayOfValues(%d)= "statore_%d,Face#%d"\n',jj-1,jj,4);
end
fprintf(fid,'%s.createBoundaryCondition(ArrayOfValues, "BoundaryCondition#2")\n',CGD);
fprintf(fid,'%s.setMagneticFluxTangential("BoundaryCondition#2")\n',CGD);
fprintf(fid,'%s.endUndoGroup()\n',CGD);

% Boundary Condition #3 - Flux Tangential - Rotor or Shaft
fprintf(fid,'REDIM ArrayOfValues(0)\n');
if (add_shaft)
    if IPM
        face_shaft = 'Face#6';
    else
        if ang_fine_r <= pi/2
            face_shaft = 'Face#3';
        else
            face_shaft = 'Face#4';
        end
    end
    fprintf(fid,'%s.selectObject("shaft,%s", infoSetSelection)\n',CGDGV,face_shaft);
    fprintf(fid,'%s.beginUndoGroup("Assign Boundary Condition")\n',CGD);
    fprintf(fid,'ArrayOfValues(0)= "shaft,%s"\n',face_shaft);
else
    fprintf(fid,'%s.selectObject("rotor,Face#6", infoSetSelection)\n',CGDGV);
    fprintf(fid,'%s.beginUndoGroup("Assign Boundary Condition")\n',CGD);
    fprintf(fid,'ArrayOfValues(0)= "rotor,Face#6"\n');
end
fprintf(fid,'%s.createBoundaryCondition(ArrayOfValues, "BoundaryCondition#3")\n',CGD);
fprintf(fid,'%s.setMagneticFluxTangential("BoundaryCondition#3")\n',CGD);
fprintf(fid,'%s.endUndoGroup()\n',CGD);
% 
% MESH DEL ROTORE
fprintf(fid,'%s.selectObject("rotor", infoSetSelection)\n',CGDGV);
fprintf(fid,'%s.beginUndoGroup("Set rotor Properties", true)\n',CGD);
if IPM
    fprintf(fid,'%s.setMaxElementSize("rotor", %d)\n',CGD,MeshElement_iron);
else
    fprintf(fid,'%s.setMaxElementSize("rotor", %d)\n',CGD,MeshElement_iron * 2);
end
fprintf(fid,'%s.endUndoGroup()\n',CGD);

% mesh dei magneti
if ~IPM
    for jj = 1:(n_magneti_polo*repetitions_r)
        fprintf(fid,'%s.selectObject("magnet_%d", infoSetSelection)\n',CGDGV,jj);
        fprintf(fid,'%s.setMaxElementSize("magnet_%d", %d)\n',CGD,jj,MeshElement_iron);
        fprintf(fid,'%s.endUndoGroup()\n',CGD);
    end
end
% %% MESH DEL TRAFERRO
% % 1 - Allinemanento Rotore - statore
% % fprintf(fid,'%s.rotateComponent (ArrayMotion, 0, 0, 0, 0, 0, 1, %d, 1)\n',CGD,(-ang_asse_avv)*180/pi);
% % 2 - mesh subdivision
if IPM
    fprintf(fid,'REDIM ArrayOfValues(0)\n');
    fprintf(fid,'ArrayOfValues(0)= "air_gap_rotor,Face#1,Edge#3"\n');
    fprintf(fid,'%s.assignMeshEdgeSubdivisions(ArrayOfValues, "Type=Uniform;Subdivisions=%d;DensityRatio=1")\n',CGD,Mac.N_airgap);
    fprintf(fid,'ArrayOfValues(0)= "air_gap_rotor,Face#2,Edge#4"\n');
    fprintf(fid,'%s.assignMeshEdgeSubdivisions(ArrayOfValues, "Type=Uniform;Subdivisions=%d;DensityRatio=1")\n',CGD,Mac.N_airgap);
else
    fprintf(fid,'REDIM ArrayOfValues(0)\n');
    % air_gap_rotor
    if ang_fine_r <= pi/2
        fprintf(fid,'ArrayOfValues(0)= "air_gap_rotor,Face#1,Edge#2"\n');
    else
        fprintf(fid,'ArrayOfValues(0)= "air_gap_rotor,Face#1,Edge#1"\n');
    end
    fprintf(fid,'%s.assignMeshEdgeSubdivisions(ArrayOfValues, "Type=Uniform;Subdivisions=%d;DensityRatio=1")\n',CGD,Mac.N_airgap/4);
    % rotor
    if ang_fine_r <= pi/2
        fprintf(fid,'ArrayOfValues(0)= "rotor,Face#2,Edge#3"\n');
    else
        fprintf(fid,'ArrayOfValues(0)= "rotor,Face#2,Edge#4"\n');
    end
    fprintf(fid,'%s.assignMeshEdgeSubdivisions(ArrayOfValues, "Type=Uniform;Subdivisions=%d;DensityRatio=1")\n',CGD,Mac.N_airgap/8);
end

% 3 - Contro Allineamento Rotore - Statore
% fprintf(fid,'%s.rotateComponent (ArrayMotion, 0, 0, 0, 0, 0, 1, %d, 1)\n',CGD,(ang_asse_avv)*180/pi);

%SAVE MACHINENAME.MN
fprintf(fid,'%s.save("%s%s.mn")',CGD,Mac.MachinePath,Mac.MachineName);

fclose(fid);
