% segmenti

segm =  [3 6             % fianco dente interno
        6 9             % raccordo squadrato fondo cava
        9 18             % fondo cava verticale
        18 12
        12 23
        23 14
        8 9
        9 10
        17 18
        4 10
        10 16];

numseg = length(segm);

% archi
arcseg=[2 3 s.AS1*acs 0.5 0         % semiapertura cava (grid 0.5deg)
        3 4 s.AS1*(1-acs) 0.5 0     % mezzo dente (grid 0.5deg)
        15 16 s.AS1 s.AS1/5 0];

numarcseg=size(arcseg,1);


