%% MG 2013/07/02 
%% Export MOGA to DXF 
% Il seguente scrip realizza l'export tra la geometria di MOGA e DXF

%%%%%%%%%%%%%%%
ROTmatr;
%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
STATmatr
%%%%%%%%%%%%%%%%
% closefemm;

%%  %%%%%
%% Salvataggio del file dxf
%%  %%%%%
mkdir(pathname_DXF);
nomefile=[pathname_DXF,filemot(1:end-4),'.dxf'];

raggi=[];
avvolgimento=[];
magneti=[];
rotore=[];

%% Creazione e salvataggio del file dxf:
DXFconv(raggi,avvolgimento,rotore2,statore,magneti,nomefile);







