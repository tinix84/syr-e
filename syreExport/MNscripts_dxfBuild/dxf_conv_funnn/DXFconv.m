% 28-05-2010 BB
% Modificati gli input per disegnare eventualmente anche l'albero del
% motore -> NO (Modifica annullata)
function DXFconv(raggi,avvolgimento,rotore,statore,magneti,nomefile,LAria)
% function DXFconv(raggi,avvolgimento,rotore,statore,magneti,nomefile,LAria,albero)

%% Scrittura del file dxf per motori
% versione 1.0
% DXFconv(rotore,statore,magneti,nomefile)
%
% RAGGI descrizione dei raggi di limitazione del disegno
% AVVOLGIMENTO descrizione dei punti dell'avvolgimento
% ROTORE    descrizione dei punti di rotore
% STATORE   descrizione dei punti di statore
% MAGNETI   descrizione dei punti dei magneti
% NOMEFILE  nome del file completo di path (system indipendet)

% Apertura file di destinazione
fid = fopen(nomefile,'w');

% Definizione sezione HEADER
header(fid)

% Definizione sezione CLASSES
classes(fid)

% Definizione sezione TABLES
tables(fid)

% Definizione sezione BLOCKS
blocks(fid)

% 28-05-2010 BB
% Modificati gli input per disegnare eventualmente anche l'albero del motore
% NO -> Modifica annullata
if nargin<7   
% Definizione sezione ENTITIES
    entities(raggi,avvolgimento,rotore,statore,magneti,fid);
else
    entities(raggi,avvolgimento,rotore,statore,magneti,fid,LAria);
end
% switch nargin
%     case 6
%        entities(raggi,avvolgimento,rotore,statore,magneti,fid);
%     case 7
%        entities(raggi,avvolgimento,rotore,statore,magneti,fid,LAria);
%     case 8
%         if ~isempty(LAria)
%           entities(raggi,avvolgimento,rotore,statore,magneti,fid,LAria,albero);
%         else
%           entities(raggi,avvolgimento,rotore,statore,magneti,fid,[],albero);
%         end
% end

% Definizione sezione OBJECTS
objects(fid)

% Chiusura del file
testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'EOF';
fprintf(fid,'%s\n',testo);

% Chiusura file di lavoro
fclose(fid);


