% disegna_rotore.m - 24 08 09 - GMP
%% DISEGNO ROTORE:
% va fatto prima di quello dello statore perchè esce temporanemaente dal
% proprio ingombro e cancellerebbe delle parti di statore

function disegna_rotore(geo,mat)

global error_code

xr = geo.xr;
Ar = geo.Ar;
p = geo.p;
g = geo.g;
pont0 = geo.pont0;
nlay = geo.nlay;
alpha = geo.alpha;
hc = geo.hc;
pont = geo.pont;
magpu = geo.magpu;
racc_pont = geo.racc_pont;
ang_pont0 = geo.ang_pont0;

magnet = mat.magnet;
steel = mat.steel;

% il rotore è gruppo 0
group = 0;

%% ELIMINO I PONT TORPPO SOTTILI
pont(pont < pont0) = 0;

% % ampiezza di un arco lungo pont0
% ang_pont0 = pont0 / xr * 180/pi;
% riquadro principale
x1 = xr; y1 = 0;
[x2,y2] = rot_point(x1,y1,180/p*pi/180);
x3 = xr-pont0; y3 = 0;
[x4,y4] = rot_point(x3,y3,180/p*pi/180);
% linea inizio polo
mi_drawline(0,0,x1,y1);
mi_selectsegment(mean([0 x1]),mean([0 y1]));
mi_setsegmentprop('APr1', 0, 1, 0, group);
% linea fine polo
mi_drawline(0,0,x2,y2);
mi_selectsegment(mean([0 x2]),mean([0 y2]));
mi_setsegmentprop('APr1', 0, 1, 0, group);
mi_clearselected;

% albero
x1 = Ar; y1 = 0;
[x2,y2] = rot_point(x1,y1,180/p*pi/180);
mi_drawarc(x1,y1,x2,y2,180/p,10);
mi_selectnode(0,0);
mi_deleteselectednodes;
mi_selectarcsegment(0,0);
mi_setarcsegmentprop(10, 'A=0', 0, 0);

% angolo barriere (viste dal centro del cerchio x0,0)
beta = 180/pi * calc_apertura_cerchio(pi/180*alpha,xr,sqrt(2)*xr);
gammam = magpu .* beta;
% raggio cerchio mediano barriera
x0 = xr * sqrt(2);
r = (x0 - xr * cos(alpha*pi/180))./(cos(beta*pi/180));

% raggi in e out delle barriere
r_all = [];
for jj = 1:nlay
    r_all = [r_all r(jj)-hc(jj)/2 r(jj)+hc(jj)/2];
end
% residuo per disegnare i ponticelli radiali
res = min(diff(sort(r_all)))/2;

% ruoto indietro - simmetria su asse x
mi_selectgroup(group); mi_moverotate(0,0,-90/p);

%% disegna le barriere di flusso (layers)
for jj = 1:2:length(r_all)
    
    % se alpha < 4 gradi non disegno il layer
    if (alpha(ceil(jj/2)) > 4.5)
        % punti su pont0 - parte interna
        [x5,y5] = rot_point(xr-pont0,0,(alpha(ceil(jj/2))+ang_pont0/2)*pi/180);
        [x6,y6] = rot_point(xr-pont0,0,(alpha(ceil(jj/2))-ang_pont0/2)*pi/180);

        %% verifiche che la barriera sia disegnabile (layer 1)
        %%%%
        % verifica 1. pont0 è centrato o no rispetto ai cerchi che delimitano il layer?
        % Succede che per alpha(1) piccolo insieme a dalpha(2) molto piccolo la
        % punta del layer si sfalsi molto rispetto alla traiettoria circolare
        % -> viene una punta che non riesce a meshare. In quel caso sposto la
        % posizione di pont0 del layer 1 più in basso. Il punto 5 bis è
        % l'intersezione di r_all(2) con la corona circolare xr - pont0.
        [x5_bis,y5_bis] = calc_intersezione_cerchi(xr-pont0,r_all(jj+1),x0);
        if x5_bis > x5
            % se il cerchio maggiore non chiappa pont0 sposto pont0 in modo da rispettare la condizione limite
            alpha_new = atan(y5_bis/x5_bis);
            [x5,y5] = rot_point(xr-pont0,0,alpha_new);
            [x6,y6] = rot_point(xr-pont0,0,alpha_new-ang_pont0*pi/180);
        end
        % verifica 2. il confine esterno di layer 1 è circolare oppure
        % verticale? Succede per alpha(1) molto piccoli, indipendentemente da
        % dalpha(2). Quando non c'è intersezione tra i cerchi tiro una
        % linea verticale.
        [x3,y3] = calc_intersezione_cerchi(xr-pont0-hc(ceil(jj/2))/2,r_all(jj),x0);
        if ~isreal(y3)
            % se non c'è intersezione (barriera piccola):
            % 1. faccio un tentativo con punta layer più piccola (1/nlay
            % invece che 1/2)
            [x3,y3] = calc_intersezione_cerchi(xr-pont0-hc(ceil(jj/2))/nlay,r_all(jj),x0);
            if ~isreal(y3)
                % 2. se ancora non funziona retta verticale invece che arco
                x3 = x6;
                % y3 = pont(ceil(jj/2))/2;
                y3 = 0;
            end
        end
        [x4,y4] = calc_intersezione_cerchi(xr-pont0-hc(ceil(jj/2))/2,r_all(jj+1),x0);
        mi_addnode(x3,y3);
        if y3 >0
            mi_addnode(x3,-y3);
        end
        mi_addnode(x4,y4);
        mi_addnode(x4,-y4);

        %% arco -hc/2
        [alpha3,r3] = cart2pol(x3,y3);
        beta3 = calc_apertura_cerchio(alpha3,r3,x0);
        mi_drawarc(x3,y3,x3,-y3,2*beta3*180/pi,5);
        %% arco +hc/2
        [alpha4,r4] = cart2pol(x4,y4);
        beta4 = calc_apertura_cerchio(alpha4,r4,x0);
        mi_drawarc(x4,y4,x4,-y4,2*beta4*180/pi,5);

        % punti su pont0 (già calcolati all'inizio)
        mi_addnode(x5,y5);
        mi_addsegment(x5,y5,x4,y4);
        mi_addnode(x5,-y5);
        mi_addsegment(x5,-y5,x4,-y4);
        mi_addnode(x6,y6);
        mi_addsegment(x6,y6,x3,y3);
        mi_addnode(x6,-y6);
        mi_addsegment(x6,-y6,x3,-y3);
        mi_addsegment(x5,y5,x6,y6);
        mi_addsegment(x5,-y5,x6,-y6);

        % magneti
        no_air = 0;
        if gammam(ceil(jj/2))>0
            % se ci sono i magneti, i punti teorici sono 7 e 8
            [x7,y7] = calc_punto_magnete(r_all(jj),gammam(ceil(jj/2))*pi/180,sqrt(2) * xr);
            [x8,y8] = calc_punto_magnete(r_all(jj+1),gammam(ceil(jj/2))*pi/180,sqrt(2) * xr);
            
            %% ECCEZIONI VERSO LA PUNTA DEL LAYER
            % 1. il magnete sconfina nella punta del layer da sotto (supera il
            % punto 3)
            if (y7>y3)
                % il magnete supera anche il punto 6 -> tutto magnete -> no_air
                [alpha7,r7] = cart2pol(x7,y7);
                beta7 = calc_apertura_cerchio(alpha7,r7,x0);
                [alpha6,r6] = cart2pol(x6,y6);
                beta6 = calc_apertura_cerchio(alpha6,r6,x0);
                if (beta7>beta6)
                    x7 = x6; y7 = y6;
                    x8 = x5; y8 = y5;
                    no_air = 1;
                else
                    retta1 = [x3 y3 x6 y6];
                    retta2 = sqrt(2)*xr*[1 0 0 tan(gammam(ceil(jj/2))*pi/180)];
                    % evito la condizione singolare retta1 verticale
                    if x3 ~= x6
                        % non verticale
                        [x7,y7] = calc_intersezione_rette(retta1,retta2);
                    else
                        % verticale
                        x7 = x3;
                        % coefficienti retta2
                        c = (retta2(4) - retta2(2))/(retta2(3)-retta2(1));
                        d = retta2(2) - c * retta2(1);
                        % intersezione con retta1 verticale
                        y7 = c*x7 + d;
                    end
                    x3 = x7; y3 = y7;
                end
            end
            % 2. il magnete sconfina nella punta del layer da sotto (supera il
            % punto 4)
            if (y8>y4)
                [alpha8,r8] = cart2pol(x8,y8);
                beta8 = calc_apertura_cerchio(alpha8,r8,x0);
                [alpha5,r5] = cart2pol(x5,y5);
                beta5 = calc_apertura_cerchio(alpha5,r5,x0);
                if (beta8>=beta5)
                    % tutto magnete
                    x7 = x6; y7 = y6;
                    x8 = x5; y8 = y5;
                    no_air = 1;
                else
                    retta1 = [x5 y5 x4 y4];
                    retta2 = sqrt(2)*xr*[1 0 0 tan(gammam(ceil(jj/2))*pi/180)];
                    [x8,y8] = calc_intersezione_rette(retta1,retta2);
                    x4 = x8; y4 = y8;
                end
            end
            % ora i punti 7 e 8 sono definiti
            % aspetto i pronticelli prima di tracciare il segmento 7-8
        end
        %% ponticelli radiali
        % punto esterno e punto interno della parte rettilinea
        x10 = sqrt(2)*xr - r_all(jj) - racc_pont;
        x11 = sqrt(2)*xr - r_all(jj+1) + racc_pont;
        % se i due punti si incrociano (ponticello piccolo rispetto al
        % raccordo) elimino il ponticello
        if x11 < x10
            % lunghezza pont (CON raccordi)
            lpont = x10 - x11;
            hpont = pont(ceil(jj/2));
            y10 = hpont/2;
            y11 = y10;
            % linea orizzontale
            mi_drawline(x11,y11,x10,y10);
            mi_drawline(x11,-y11,x10,-y10);
            % raccordo esterno: verifico che non cozzi con i magneti
            x12 = x10 + racc_pont; y12 = y10 + racc_pont;
            if y12 > y7
                % il magnete x è troppo piccolo, non lo disegno
                y7 = 0; y8 = 0;
            end
            mi_drawline(x10,y10,x12+res,y12+res);
            mi_drawline(x10,-y10,x12+res,-y12-res);
            % cancello i nodi di troppo
            mi_selectnode(x12+res,y12+res);
            mi_selectnode(x12+res,-y12-res);
            mi_deleteselectednodes;

            % raccordo interno
            x13 = x11 - 1.5 * racc_pont; y13 = y10+1.5*racc_pont;
            mi_drawline(x11,y11,x13,y13);
            mi_drawline(x11,-y11,x13,-y13);
             
            % cancello i nodi di troppo
            mi_selectnode(x13,y13);
            mi_selectnode(x13,-y13);
            mi_deleteselectednodes
            % cancello gli archetti di confine
            mi_selectarcsegment(x11,0);
            mi_deleteselectedarcsegments;
            if x3 ~= x6
                mi_selectarcsegment(x10,0);
                mi_deleteselectedarcsegments;
            else
                mi_selectnode(x3,0);
                mi_deleteselectednodes;
            end
        else
            % non faccio il ponticello disegno solo una linea di separazione
            mi_drawline(x11- racc_pont,0,x3,0);
            hpont = 0;
        end
        % assegno aria
        if no_air == 0
            % assegna aria - nel baricentro della punta, delimitata dai nodi 3-4-5-6
            x = mean([x3 x4 x5 x6]); y = mean([y3 y4 y5 y6]);
            mi_addblocklabel(x,y); mi_addblocklabel(x,-y);
            mi_selectlabel(x,y); mi_selectlabel(x,-y);
            mi_setblockprop('Air', 1, 0, 'None', 0, 0, 1);
            mi_clearselected;
        end
       
        % OK magneti:
        if (gammam(ceil(jj/2))>0 && (y7 >0))
            % - disegno segmento 7-8
            mi_addnode(x7,y7); mi_addnode(x7,-y7);
            mi_addnode(x8,y8); mi_addnode(x8,-y8);
            mi_addsegment(x7,y7,x8,y8);
            mi_addsegment(x7,-y7,x8,-y8);
            % - assegno materiale
            x = sqrt(2) * xr - r(ceil(jj/2));
            y = hpont;
            % eccezione per magneti molto piccoli
            if y > y7
                y = mean([hpont/2 y7]);
            end
            % eccezione quando non c'è pont
            if hpont == 0
                y = y7/30;
            end
            mi_addblocklabel(x,y); mi_selectlabel(x,y);
            mi_setblockprop(magnet, 1, 0, 'None', -gammam(ceil(jj/2))/2 + 180, 0, 1);
            mi_clearselected;
            mi_addblocklabel(x,-y); mi_selectlabel(x,-y);
            mi_setblockprop(magnet, 1, 0, 'None', gammam(ceil(jj/2))/2 + 180, 0, 1);
            mi_clearselected;
        end
    end
end
mi_clearselected;

% disegno arco esterno
% singolo arco affacciato al traferro
x = xr * cos(90/p * pi/180);
y = xr * sin(90/p * pi/180);
mi_drawarc(x,-y,x,y,180/p,1);
% % assegno magneti e ferro
% for jj = 1:nlay
%     if gammam(jj)>0
%         keyboard
%         x = sqrt(2) * xr - r(jj);
%         y = pont(jj);
%         % eccezione per magneti molto piccoli
%         if y > y7
%             y = mean([pont(jj)/2 y7]);
%         end
%         mi_addblocklabel(x,y); mi_selectlabel(x,y);
%         mi_setblockprop(magnet, 1, 0, 'None', -gammam(jj)/2 + 180, 0, 1);
%         mi_clearselected;
%         mi_addblocklabel(x,-y); mi_selectlabel(x,-y);
%         mi_setblockprop(magnet, 1, 0, 'None', gammam(jj)/2 + 180, 0, 1);
%         mi_clearselected;
%     end
% end
% mi_clearselected;
% assegno ferro
mi_addblocklabel(xr - pont0,0);
mi_selectlabel(xr - pont0,0);
mi_setblockprop(steel, 1, 0, 'None', 0, 0, 1);

% ruoto avanti
mi_selectgroup(group); mi_moverotate(0,0,90/p);

% gap rotore
x1 = xr; y1 = 0;
x2 = xr+1/3*g; y2 = 0;
mi_drawline(x1,y1,x2,y2);
mi_selectsegment(mean([x1 x2]),mean([y1 y2]));
mi_setsegmentprop('APg1', 0, 1, 0, group);

[x1,y1] = rot_point(x1,y1,180/p*pi/180);
[x2,y2] = rot_point(x2,y2,180/p*pi/180);
mi_drawline(x1,y1,x2,y2);
mi_selectsegment(mean([x1 x2]),mean([y1 y2]));
mi_setsegmentprop('APg1', 0, 1, 0, group);

[x1,y1] = rot_point(xr+1/6*g,0,90/p*pi/180);
mi_addblocklabel(x1,y1);
mi_selectlabel(x1,y1);
mi_setblockprop('Air', 1, 0, 'None', 0, group, 1);

% mi_saveas(filename) % saves the file with name ’filename’.
% mi_zoomnatural
% pause




