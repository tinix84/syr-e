% disegna_motore.m - 24 08 09 - GMP
% disegna il motore in posizione zero

% input:
% - geo : geometria
% - mat : materiali

% output:
% - mot0.fem (th_m = 0, i123 = 0)

%% modifica 10 gennaio 2010 - gp
% - porto fuori geo anche da disegna_statore() per salvare i punti x_fe,
% y_fe di statore in geo_mot_temp

error_code = 0;

% carica il file vuoto (controllare che abbia tutti i materiali caricati)
%% COMMENTO PARFOR
%    h_temp=actxserver('femm.ActiveFEMM')
%    callfemm([ 'setcurrentdirectory(' , quote(pwd) , ')'],h_temp)  
 filename = 'mot0.fem';
% opendocument('empty_case.fem',h_temp);
 opendocument('empty_case.fem');

mi_probdef(0,'millimeters','planar',1e-8,geo.l,30);

% fasatura meccanica
th_m0 = 7.5;
geo.th_m0=th_m0;
% boundary
% A = 0 sulle superfici interna ed esterna
mi_addboundprop('A=0', 0, 0, 0, 0, 0, 0, 0, 0, 0);

% anti-periodicità (2 x rotore + 2 x statore + 3 x airgap + AP move)
mi_addboundprop('APr1', 0, 0, 0, 0, 0, 0, 0, 0, 5);
mi_addboundprop('APr2', 0, 0, 0, 0, 0, 0, 0, 0, 5);

mi_addboundprop('APg1', 0, 0, 0, 0, 0, 0, 0, 0, 5);
mi_addboundprop('APg2', 0, 0, 0, 0, 0, 0, 0, 0, 5);
mi_addboundprop('APg3', 0, 0, 0, 0, 0, 0, 0, 0, 5);
mi_addboundprop('APmove', 0, 0, 0, 0, 0, 0, 0, 0, 5);

mi_addboundprop('APs1', 0, 0, 0, 0, 0, 0, 0, 0, 5);
mi_addboundprop('APs2', 0, 0, 0, 0, 0, 0, 0, 0, 5);


%keyboard
geo = disegna_rotore_semicerchi(geo,mat);          % NB: Il disegno del rotore va fatto prima di quello dello statore perchè esce temporanemaente dal proprio 
%                                                    % ingombro e
% %keyboard                                                  % cancellerebbe delle
%                                                    % parti di statore
% 
% % posizione angolare (serve per la fasatura con l'avvolgimento)
mi_selectgroup(2); 
mi_moverotate(0,0,th_m0);

% disegna statore e salva mot0
geo = disegna_statore(geo,mat);

% disegna bordo mobile (group 20)
disegna_bordo_mobile(geo,th_m0);
keyboard
mi_saveas(filename); % saves the file with name ’filename’.
% closefemm;




