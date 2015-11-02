%% VALUTAZIONE DELLA COPPIA PER (nsim-1) DIVERSE POSIZIONI ANGOLARI DEL ROTORE (COMPRESE TRA 0 E XDEG) _ FASE DELLA CORRENTE COSTANTE (GAMMA)

% Input: - nsim; - id, iq;
% Output: - vettore soluzione SOL (nsim x 6); - Temp\sim_gamma_numerosimulazione.fem; - Temp\sim_gamma.mat (memorizza SOL)

% Legenda delle colonne di SOL -> 1 - gamma,  2:3 - id, iq,  4:5 - fd, fq,
% 6 - coppia 7-end Bx By in vari punti (Pfe)

%% modifica 10 gennaio 2010 - gp
% - aggiunto tipo_valutazione per distinguere tra: gambjfitness, eval_x,
% valuta_singola_macchina

function SOL = simula_xdeg(geo,nsim,xdeg,io,gamma)

global tipo_valutazione

th0=geo.th_0;
p = geo.p;
xr = geo.xr;
gap = geo.g;
ns = geo.ns;
pc = 360/(ns*p)/2;
nlay = geo.nlay;
r_fe = geo.r_fe;
hf = geo.hf;
beta_f = geo.beta_f;
l = geo.l;
x_fe = geo.x_fe;
y_fe = geo.y_fe;
id = -io * sin(gamma * pi/180);
iq = io * cos(gamma * pi/180);
Hc = geo.Hc;
% keyboard
SOL = [];

% th = linspace(0,xdeg,nsim);
% la seguente linea introduce una variabilità delle posizioni simulate in
% modo da evitare di prendere tutti i passaggi per lo zero della terza
% armonica. Si aggiunge un numero casuale tra zero e 2.5 che corrisponde ad
% un quarto di periodo della terza armonica
% th = linspace(0,xdeg,nsim)+xdeg/12*rand;

% modifiche del 9 gennaio
% il passo di simulazione dovrevve essere xdeg/(nsim-1)
% il +0.5 diminuisce la varianza di risultati in presenza di armoniche
% e permette di beccare la terza anche con tre punti
if tipo_valutazione == 'gambj'
    passo=xdeg/(nsim-1+0.5);
    % la variabilità del punto iniziale è di mezzo passo o di un passo
    caso=1*passo*rand; % un passo 1* o mezzo passo 0.5*
else
    passo=xdeg/(nsim-1);
    caso=0;
end
teta=[0:passo:xdeg]+caso;

% prendo solo i primi nsim-1 punti
th=th0+[teta(1:nsim-1) teta(1)];
% keyboard
% Calcolo subito le correnti di fase per le nsim simulazioni
% le prime nsim-1 simulazioni usano gamma e ruotano il rotore
% la nsim-esima simulazione è con gamma = 90°
% per ora in totale 6 simulazioni
%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%
%%%% AGGIUSTARE PER IL CASO DI TANTI GAMMA %%%%%%%%%%%%
%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%
for ij=1:nsim
    if ij==nsim
        th_m = th(ij)/p;
        i123 = dq2abc(-io,0,th(ij)*pi/180);
        i1_tmp(ij) = i123(1); i2_tmp(ij) = i123(2); i3_tmp(ij) = i123(3);
    else
        th_m = th(ij)/p;
        i123 = dq2abc(id,iq,th(ij)*pi/180);
        i1_tmp(ij) = i123(1); i2_tmp(ij) = i123(2); i3_tmp(ij) = i123(3);
    end
end
% keyboard
%% sceglie se è in matlabpool o no
isOpen = matlabpool('size');
if isOpen > 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% ciclo PARFOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % tic %decommentare per verificare durata calcoli
    
    names_o{1} = 'run1';names_o{2} = 'run2';names_o{3} = 'run3';names_o{4} = 'run4';
    names_o{5} = 'run5';names_o{6} = 'run6';names_o{7} = 'run7';names_o{8} = 'run8';
    names_o{9} = 'run9';names_o{10} = 'run10';names_o{11} = 'run11';names_o{12} = 'run12';
    names_o{13} = 'run13';names_o{14} = 'run14';names_o{15} = 'run15';names_o{16} = 'run16';
    names_o{17} = 'run17';names_o{18} = 'run18';names_o{19} = 'run19';names_o{20} = 'run20';
    
    parfor jj = 1:nsim %% 9.10.11 aggiungo una simulazione per incorporare anche sol90
        
        %tmp_file=tempname;
        tmp_fem=[names_o{jj} '.fem'];
        copyfile('mot0.fem',tmp_fem);
        h_temp=actxserver('femm.ActiveFEMM');
        callfemm_parfor([ 'setcurrentdirectory(' , quote(pwd) , ')'],h_temp);
        opendocument_parfor(tmp_fem,h_temp);
        
        th_m = (th(jj) - th0)/p;
        %     % correnti di fase
        i1 = i1_tmp(jj);i2 = i2_tmp(jj);i3 = i3_tmp(jj);
        %
        %     % nome caso
        %     filename = ['Temp\sim_' num2str(io) '_' num2str(gamma) '_n' num2str(th_m)];
        %
        %     %opendocument('mot0.fem');
        % assegna le correnti
        mi_modifycircprop('fase1',1,i1);
        mi_modifycircprop('fase1n',1,-i1);
        mi_modifycircprop('fase2',1,i2);
        mi_modifycircprop('fase2n',1,-i2);
        mi_modifycircprop('fase3',1,i3);
        mi_modifycircprop('fase3n',1,-i3);
        % assegna il magnete
        mi_modifymaterial('Bonded-Magnet',3,Hc);
        % cancello bordo mobile
        mi_selectgroup(20), mi_deleteselectedarcsegments;
        % posiziono il rotore
        mi_selectgroup(2), mi_moverotate(0,0,th_m);
        
        % disegna bordo mobile (group 20)
        disegna_bordo_mobile(geo,th_m);
        
        %     mi_saveas('mot_temp.fem');
        
        mi_analyze(1);
        
        % scarica i risultati
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%% INCOLLO QUI POST_PROC %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% FC 9.10.11 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %     post_proc_parfor;
        
        %% modifica 10 gennaio 2010 - gp
        % - scarico Bx By anche nel dente di statore
        % - organizzazione dei vettori x_fe, y_fe
        % x_fe(1:nlay) : barriere di rotore, vicino alle punte tonde
        % x_fe(nlay+1:nlay+2) : espansione dente statore (speculari)
        % x_fe(nlay+3:nlay+4) : dente statore (speculari)
        % x_fe(nlay+5) : giogo statore
        
        mi_loadsolution_parfor(h_temp);
        
        % evaluate flux
        temp_out = mo_getcircuitproperties('fase1');
        temp_out = temp_out - mo_getcircuitproperties('fase1n');
        f1 = temp_out(3) * 2 * p;
        temp_out = mo_getcircuitproperties('fase2');
        temp_out = temp_out - mo_getcircuitproperties('fase2n');
        f2 = temp_out(3) * 2 * p;
        temp_out = mo_getcircuitproperties('fase3');
        temp_out = temp_out - mo_getcircuitproperties('fase3n');
        f3 = temp_out(3) * 2 * p;
        
        % evaluate torque
        % T1 - linea vicina al rotore
        x = xr + gap*1/6;
        ang0 = th_m; ang1 = 180/p + th_m;
        [x1,y1] = rot_point(x,0,ang0*pi/180);
        [x2,y2] = rot_point(x,0,ang1*pi/180);
        mo_addcontour(x1,y1);
        mo_addcontour(x2,y2);
        mo_bendcontour(180/p,0.5);
        
        T1 = mo_lineintegral(4);
        T1 = T1(1) * 2 * p;
        mo_clearcontour();
        
        % T2 - linea vicina allo statore
        x = xr + gap*5/6;
        ang0 = -pc; ang1 = 180/p-pc;
        
        [x1,y1] = rot_point(x,0,ang0*pi/180);
        [x2,y2] = rot_point(x,0,ang1*pi/180);
        mo_addcontour(x1,y1);
        mo_addcontour(x2,y2);
        mo_bendcontour(180/p,0.5);
        T2 = mo_lineintegral(4);
        T2 = T2(1) * 2 * p;
        mo_clearcontour();
        
        % T3 - linea mediana
        x = xr + gap*1/2;
        ang0 = -pc; ang1 = 180/p-pc;
        
        [x1,y1] = rot_point(x,0,ang0*pi/180);
        [x2,y2] = rot_point(x,0,ang1*pi/180);
        mo_addcontour(x1,y1);
        mo_addcontour(x2,y2);
        mo_bendcontour(180/p,0.5);
        T3 = mo_lineintegral(4);
        T3 = T3(1) * 2 * p;
        mo_clearcontour();
        
        %% PERDITE NEL FERRO
        % Pfe ROT - valuto il campo nelle guide di flusso
        % - la guida più superficiale (numero 1) non viene considerata
        % - la guida più interna (vicino all'albero) è trattata in modo particolare (punto 4 in disegna_rotore.m) perchè
        % non ha spessore uniforme (nella zona dell'albero il campo è molto debole e non significativo)
        
%         % Punto 0: centro guida sull'asse di simmetria del polo
%         x_fe0 = geo.x0 - r_fe;
%         ang_temp = (90/p + th_m)*pi/180;
%         [x_fe0,y_fe0] = rot_point(x_fe0,0,ang_temp);
%         
%         % Bxy rotore
%         temp_Bx0 = zeros(1,nlay); temp_By0 = temp_Bx0;
%         temp_Bx1 = temp_Bx0; temp_By1 = temp_Bx0;
%         temp_Bx2 = temp_Bx0; temp_By2 = temp_Bx0;
%         % Bxy statore
%         temp_Bx3 = zeros(1,ns/2); temp_By3 = temp_Bx3;
%         
%         % Punto 1: centro guida vicino al traferro
%         [x_fe1,y_fe1] = rot_point(x_fe(1:nlay),y_fe(1:nlay),ang_temp);
%         
%         % Punto 2: specchio di punto 1
%         [x_fe2,y_fe2] = rot_point(x_fe(1:nlay),-y_fe(1:nlay),ang_temp);
%         
%         % scarico il campo di rotore
%         for ii = 1:nlay
%             temp = mo_getpointvalues(x_fe0(ii),y_fe0(ii));
%             temp_Bx0(ii) = temp(2);
%             temp_By0(ii) = temp(3);
%             temp = mo_getpointvalues(x_fe1(ii),y_fe1(ii));
%             temp_Bx1(ii) = temp(2);
%             temp_By1(ii) = temp(3);
%             temp = mo_getpointvalues(x_fe2(ii),y_fe2(ii));
%             temp_Bx2(ii) = temp(2);
%             temp_By2(ii) = temp(3);
%         end
%         [Bx0,By0] = rot_point(temp_Bx0,temp_By0,-ang_temp);
%         [Bx1,By1] = rot_point(temp_Bx1,temp_By1,-ang_temp);
%         [Bx2,By2] = rot_point(temp_Bx2,temp_By2,-ang_temp);
%         
%         % Pfe STAT - valuto il campo nel primo dente
%         
%         % valuto tutti i denti di statore (ns/2)
%         % x_fe3 = zeros(1,ns/2); y_fe3 = x_fe3;
%         x_fe3 = x_fe(nlay+3);
%         y_fe3 = y_fe(nlay+3);
%         
%         for ii = 1:ns/2
%             ang_cava = pc * pi/180 * (ii-1);
%             [temp_x,temp_y] = rot_point(x_fe3,y_fe3,ang_cava);
%             temp = mo_getpointvalues(temp_x,temp_y);
%             temp_Bx3(ii) = temp(2);
%             temp_By3(ii) = temp(3);
%         end
%         Bx3 = temp_Bx3;
%         By3 = temp_By3;
        
        mo_close, mi_close
        
        Bx0=0; By0=0; Bx1=0; By1=0; Bx2=0; By2=0; Bx3=0; By3=0;
        % flussi dq
        fdq = abc2dq(f1,f2,f3,th(jj)*pi/180);
        
        sol = [th(jj) id iq fdq(1) fdq(2) mean([T1,T2,T3]) Bx0 By0 Bx1 By1 Bx2 By2 Bx3 By3]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% FINE DI POST PROC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        SOL_PARFOR_temp(jj,:) = sol;
        %     ripple_pu = abs(std(SOL_PARFOR(:,end))/mean(SOL_PARFOR(:,end)));
        %    mi_loadsolution_parfor(h_temp);
        %     pause(1)
        h_temp.delete;
        %   copyfile(tmp_fem,[names_o{jj} '.fem']);
        delete(tmp_fem);
        ans_temp=strrep(tmp_fem,'.fem','.ans');
        %    copyfile(ans_temp,[names_o{jj} '.ans']);
        delete(ans_temp);
    end
    
    % toc %decommentare per verificare durata calcoli
    
    % ripple_pu = abs(std(SOL_PARFOR_temp(:,6))/mean(SOL_PARFOR_temp(:,6)));
    
    SOL = SOL_PARFOR_temp;
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% ciclo for %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SOL=[]; %riazzero SOL

    for jj = 1:nsim
        th_m = (th(jj) - th0)/p;
        % correnti di fase
        i1 = i1_tmp(jj);i2 = i2_tmp(jj);i3 = i3_tmp(jj);
        %%     i123 = dq2abc(id,iq,th(jj)*pi/180);
        %%     i1 = i123(1); i2 = i123(2); i3 = i123(3);
        
        % nome caso
        filename = ['Temp\sim_' num2str(io) '_' num2str(gamma) '_n' num2str(th_m)];
        
        opendocument('mot0.fem');
        % assegna le correnti
        mi_modifycircprop('fase1',1,i1);
        mi_modifycircprop('fase1n',1,-i1);
        mi_modifycircprop('fase2',1,i2);
        mi_modifycircprop('fase2n',1,-i2);
        mi_modifycircprop('fase3',1,i3);
        mi_modifycircprop('fase3n',1,-i3);
        % assegna il magnete
        mi_modifymaterial('Bonded-Magnet',3,Hc);
        % cancello bordo mobile
        mi_selectgroup(20), mi_deleteselectedarcsegments;
        % posiziono il rotore
        mi_selectgroup(2), mi_moverotate(0,0,th_m);
        % disegna bordo mobile (group 20)
        disegna_bordo_mobile(geo,th_m);
        mi_saveas('mot_temp.fem');
        mi_analyze(1);
        % scarica i risultati
%         keyboard
        post_proc;
        SOL = [SOL; sol];
    end
    % ripple_pu = abs(std(SOL(:,6))/mean(SOL(:,6)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       FINE VECCHIO CICLO FOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%% VERIFICA CHE I DUE CICLI FOR E PARFOR DIANO LO STESSO RISULTATO %%
% ERROR = SOL-SOL_PARFOR_temp
