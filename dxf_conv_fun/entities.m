% 28-05-2010 BB
% Modificati gli input per disegnare eventualmente anche l'albero del
% motore -> NO Modifica annullata
function entities(raggi,avvolgimento,rotore,statore,magneti,fid,LAria)
% function entities(raggi,avvolgimento,rotore,statore,magneti,fid,LAria,albero)

global PC

%% Scrittura dei punti che caratterizzano i diversi layer
%color 0 = nero
%color 1 = rosso
%color 2 = giallo
%color 3 = verde
%color 4 = ciano
%color 5 = blu
%color 6 = magenta

% 22-06-2011 modifica dei file dei colori per poter cambiare la stampa del
% DXF
% Recupero dei colori
% C_sta=PC.DXF.Color.sta;
% C_rot=PC.DXF.Color.rot;
% C_mag=PC.DXF.Color.mag;
% C_avv=PC.DXF.Color.avv;
% C_rag=PC.DXF.Color.rag;
% C_ari=PC.DXF.Color.ari;
%% 2013/06/28 MG modifica provvisoria (ripristinare in seguito le righe precedenti)
C_sta='0';
C_rot='0';
C_mag='0';
C_avv='0';
C_rag='0';
C_ari='0';

% generazione sezione ENTITIES
% Scrittura intestazione ENTITIES
testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'SECTION';
fprintf(fid,'%s\n',testo);
testo = '2';
fprintf(fid,'%s\n',testo);
testo = 'ENTITIES';
fprintf(fid,'%s\n',testo);

% Riconoscimento della matrice (statore)
if ~isempty(statore)
    dimobj = size(statore);
    righe = dimobj(1);
    layer = '1';
    color = C_sta;
    TipoLinea = 'Continuous';
    
    for k = 1:righe
        if statore(k,7) == 0
            linea(statore,layer,color,TipoLinea,k,fid);
        elseif abs(statore(k,7)) == 1
            arco(statore,layer,color,TipoLinea,k,fid);
        else
            break
        end
    end
end

% Riconoscimento della matrice (rotore)
if ~isempty(rotore)
    dimobj = size(rotore);
    righe = dimobj(1);
    layer = '2';
    color = C_rot;
    TipoLinea = 'Continuous';
    
    for k = 1:righe
        if rotore(k,7) == 0
            linea(rotore,layer,color,TipoLinea,k,fid);
        elseif abs(rotore(k,7)) == 1
            arco(rotore,layer,color,TipoLinea,k,fid);
        else
            break
        end
    end
end

% Riconoscimento della matrice (magneti)
if ~isempty(magneti)
    dimobj = size(magneti);
    righe = dimobj(1);
    layer = '2';
    color = C_mag;
    TipoLinea = 'Continuous';
    
    for k = 1:righe
        if magneti(k,7) == 0
            linea(magneti,layer,color,TipoLinea,k,fid);
        elseif abs(magneti(k,7)) == 1
            arco(magneti,layer,color,TipoLinea,k,fid);
        else
            break
        end
    end
end


%riconoscimento della matrice (Rame)
if ~isempty(avvolgimento)
    dimobj = size(avvolgimento);
    righe = dimobj(1);
    layer = '3';
    color = C_avv;
    TipoLinea = 'Continuous';
    
    for k = 1:righe
        if avvolgimento(k,7) == 0
            linea(avvolgimento,layer,color,TipoLinea,k,fid);
        elseif abs(avvolgimento(k,7)) == 1
            arco(avvolgimento,layer,color,TipoLinea,k,fid);
        else
            break
        end
    end
end

%riconoscimento della matrice (Raggi)
if ~isempty(raggi)
    dimobj = size(raggi);
    righe = dimobj(1);
    layer = '4';
    color = C_rag;
    TipoLinea = 'Continuous';
    
    for k = 1:righe
        if raggi(k,7) == 0
            linea(raggi,layer,color,TipoLinea,k,fid);
        elseif abs(raggi(k,7)) == 1
            arco(raggi,layer,color,TipoLinea,k,fid);
        else
            break
        end
    end
end

% 28-05-2010 BB
% Gli input da raggi (il primo) a fid (l'ultimo prima degli eventuali
% LAria e albero) sono 6. Quindi devo controllare se le matrici LAria e albero
% sono oppure no vuote se gli input sono pi� di sei

% NO -> Modifica annullata

if nargin>7
% if nargin>6
    %riconoscimento della matrice (LAria)
    if ~isempty(LAria)
        dimobj = size(LAria);
        righe = dimobj(1);
        layer = '5';
        color = C_ari;
        TipoLinea = 'Continuous';
        
        for k = 1:righe
            if LAria(k,7) == 0
                linea(LAria,layer,color,TipoLinea,k,fid);
            elseif abs(raggi(k,7)) == 1
                arco(LAria,layer,color,TipoLinea,k,fid);
            else
                break
            end
        end
    end
    % 28-05-2010 BB
    % Modificati gli input per disegnare eventualmente anche l'albero del
    % motore
    %riconoscimento della matrice (albero)
    % if nargin==8 && (~isempty(albero))
    %     dimobj = size(albero);
    %     righe = dimobj(1);
    %     layer = '6';
    %     color = '2';
    %     TipoLinea = 'Continuous';
        
    %     for k = 1:righe
    %         if albero(k,7) == 0
    %             linea(albero,layer,color,TipoLinea,k,fid);
    %         elseif abs(albero(k,7)) == 1
    %             arco(albero,layer,color,TipoLinea,k,fid);
    %         else
    %             break
    %         end
    %     end
    % end
end

testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'ENDSEC';
fprintf(fid,'%s\n',testo);


