%% VERSIONE 20 11 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script utilizzato in BuildNewMachine.m

% Input:    Parametri del motore (ParMachine.m)
%           Matrici dei nodi ('motore'.mat)
% Output:   statore.vbs (script x la costruzione in MN dello statore )
%           rotore.vbs  (script x la costruzione in MN del rotore )
%           'motore'.mn (modello del motore in magnet)



% stringhe di uso comune
CGD='Call getDocument()';
CGDGV='Call getDocument().getView()';
IMCUS='infoMakeComponentUnionSurfaces';
IMCRV='infoMakeComponentRemoveVertices';
ITIS=' infoToggleInSelection';
IATS='infoAddToSelection';

% VBS script 1
build_statore_vbs_new;
% VBS script 2
build_rotore_vbs_new;
