%%

switch geo.RotType
    
    case 'Circular'
        
        if ~exist('X3')
            X3 = xxD2k;
            Y3 = yyD2k;
            X4 = xxD1k;
            Y4 = yyD1k;
        end
        
        A = [];
        rG = [];  % Distanza da (0,0) dei baricentri dei volumi appesi ai ponticelli
        tmp_bary = [];
        for j = 1 : nlay
            if j == 1
                [x,y] = calc_intersezione_cerchi(r,r_all(1),x0);
                dist = x0 - x; % distanza da x0 al punto di intersezione tra le circonferenze di raggio r e raggio r1 e centro x0
                theta = atan2(y,dist);
                theta2 = atan2(y,x);
                A1 = 0.5*(r^2*theta2 - abs(0.5*det([x y 1; 0 0 1; x -y 1])));
                A2 = 0.5*(r_all(j)^2*theta - abs(0.5*det([x y 1; x0 0 1; x -y 1])));
                A(j) = A1+A2;
                bary1 = (2*r*(sin(theta2))^3)/(3*(theta2-sin(theta2)*cos(theta2)));
                bary2 = x0-(2*r_all(1)*(sin(theta))^3)/(3*(theta-sin(theta)*cos(theta)));
                tmp_bary(j,:) = (2*A1*bary1 + 2*A2*bary2)/(2*(A1+A2));
                
            else
                [x,y] = calc_intersezione_cerchi(r,r_all(2*j-1),x0);
                dist = x0 - x;
                theta = atan2(y,dist);
                A(j) = (r_all(2*j-1)^2 - r_all(2*j-2)^2)*theta*0.5;
                theta2 = atan2(y,x);
                tmp_bary(j,:) = x0 - (2/3*sin(theta)/theta*(r_all(2*j-1)^3 - r_all(2*j-2)^3)/(r_all(2*j-1)^2 - r_all(2*j-2)^2));
                
                
            end
        end
        Afe = cumsum(A);
        %MarcoP
        rG(1) = tmp_bary(1,1);     % baricentro della regione di ferro sorretta dal ponticello a I
        for ii = 2 : nlay
            rG(ii)=(Afe(ii-1)*rG(ii-1)+A(ii)*tmp_bary(ii,1))/Afe(ii);  % baricentri delle regioni di ferro sorrette dalle due U
        end
        
        M_Fe = 2*Afe*l * 1e-9 * rhoFE ;   % massa ferro appeso ai ponticelli
        
        F_centrifuga = M_Fe .* rG/1000 *  (nmax * pi/30)^2;
        
        if geo.radial_ribs_eval == 0
            pont = F_centrifuga/(sigma_max * l);    % mm
        else
            pont = geo.pont;
        end
        
        for jj = 1 : nlay
            if (pont(jj) < pont0) % non disegno i ponticelli radiali il cui spessore è minore della tolleranza di lavorazione per gli altri tipi di rotore
                pont(jj) = 0;
            end
        end
        
        %%% DISEGNO PINTICELLI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for jj = 1 : 2 : length(r_all)
            % DISEGNO DEI PONTICELLI RADIALI
            % keyboard
            x10 = x0 - r_all(jj) - racc_pont;                           % Coordinata x Nodo 10 ((Ponticello radiale_vertice in alto a dx)=(Punto
            % in basso a dx del raccordo))
            x11 = x0 - r_all(jj+1) + racc_pont;                         % Coordinata x Nodo 11 ((Ponticello radiale_vertice in alto a sx)=(Punto
            % in basso a sx del raccordo))
            if (x11 < x10) && (pont(ceil(jj/2))>0)                      % Se i due punti non si incrociano (raccordo non troppo grande rispetto alla
                % larghezza della barriera), e se lo spessore estratto dall'algoritmo per il
                % ponticello non è troppo piccolo, allora procedo al disegno del ponticello
                %                 lpont = x10 - x11;                                      % lpont=lunghezza ponticello=larghezza barriera-2*racc_pont
                %                 hpont = pont(ceil(jj/2));
                %                 y10 = hpont/2;                                          % Coordinata y Nodo 10
                %                 y11 = y10;                                              % Coordinata y Nodo 11
                %                 XpontRadDx(ceil(jj/2)) = x10; YpontRadDx(ceil(jj/2)) = y10;
                %                 XpontRadSx(ceil(jj/2)) = x11; YpontRadSx(ceil(jj/2)) = y10;
                %                 x12 = x10 + 1.5 * racc_pont; y12 = y10 + 1.5*racc_pont;   % Punto 12 (serve solo per il disegno del raccordo: verrà cancellato)
                %                 XpontRadBarDx(ceil(jj/2)) = XBanqdx(ceil(jj/2)); YpontRadBarDx(ceil(jj/2)) = y10 + (XpontRadBarDx(ceil(jj/2)) - x10);
                %                 x13 = x11 - 1.5 * racc_pont; y13 = y10 + 1.5*racc_pont;   % Punto 13 (serve solo per il disegno del raccordo: verrà cancellato)
                %                 XpontRadBarSx(ceil(jj/2)) = XBanqsx(ceil(jj/2)); YpontRadBarSx(ceil(jj/2)) = y10 + (x11-XpontRadBarSx(ceil(jj/2)));
                
                hpont = pont(ceil(jj/2));
                y10 = hpont/2;                                          % Coordinata y Nodo 10
                YpontRadBarDx(ceil(jj/2)) = y10 + racc_pont;             % Coordinata y pontRadBarDx
                [x_temp,~] = intersezione_retta_circonferenza(x0,0,r_all(jj),0,YpontRadBarDx(ceil(jj/2)));
                XpontRadBarDx(ceil(jj/2)) = 2*x0 - x_temp;
                YpontRadBarSx(ceil(jj/2)) = y10 + racc_pont;             % Coordinata y pontRadBarDx
                [x_temp,~] = intersezione_retta_circonferenza(x0,0,r_all(jj+1),0,YpontRadBarSx(ceil(jj/2)));
                XpontRadBarSx(ceil(jj/2)) = 2*x0 - x_temp;
                XpontRadDx(ceil(jj/2))  = XpontRadBarDx(ceil(jj/2))  - racc_pont;
                YpontRadDx(ceil(jj/2))  = y10;
                XpontRadSx(ceil(jj/2))  = XpontRadBarSx(ceil(jj/2))  + racc_pont;
                YpontRadSx(ceil(jj/2))  = y10;
                
            else                                                         % Se invece i punti 10 e 11 si incrociano (ovvero raccordo troppo grande rispetto)
                % alla larghezza barriera), non faccio il ponticello e disegno solo una linea di
                hpont = 0;
                XpontRadBarSx(ceil(jj/2)) = XBanqsx(ceil(jj/2));
                YpontRadBarSx(ceil(jj/2)) = 0;
                XpontRadBarDx(ceil(jj/2)) = XBanqdx(ceil(jj/2));
                YpontRadBarDx(ceil(jj/2)) = 0;
                XpontRadDx(ceil(jj/2)) = NaN;
                YpontRadDx(ceil(jj/2)) = 0;
                XpontRadSx(ceil(jj/2)) = NaN;
                YpontRadSx(ceil(jj/2)) = 0;
            end
            if ceil(jj/2) == 1;
                x = x0 - rbeta(ceil(jj/2));
                y = hpont(1) * 0.6;
                if y == 0
                    y = pont0;
                end
            else
                [xtemp,ytemp] = rot_point(-rbeta(ceil(jj/2)),0,-0.5*beta(ceil(jj/2))*pi/180);
                x = x0 + xtemp;
                y = ytemp;
            end
        end
        %% PARTE AGGIUNTIVA PM  ===================================================
        r_semi = [];
        for iflag = 1 : nlay
            if XpontRadDx(iflag) > 0
                flag_PM(iflag) = 1;
            else
                flag_PM(iflag) = 0;
            end
        end
        % Area totale delle barriere di flusso
        for i = 1 : 2 : length(r_all)
            idx = ceil(i/2);
            dist = x0 - X3(idx);
            theta(idx) = atan2(Y3(idx),dist);
            r_semi(idx) = sqrt((X4(idx)-X3(idx))^2+(Y4(idx)-Y3(idx))^2)*0.5;
            if rem(i,2) ~= 0
                if flag_PM(idx) == 0
                    AreaBarr(idx) = pi*r_semi(idx)^2 + (r_all(2*idx)^2 - r_all(2*idx - 1)^2)*theta(idx);
                else
                    AreaBarr(idx) = pi*r_semi(idx)^2 + (r_all(2*idx)^2 - r_all(2*idx - 1)^2)*theta(idx) - flag_PM(idx)*(((XpontRadDx(idx)-XpontRadSx(idx))*2*YpontRadDx(idx)) + (YpontRadDx(idx) + YpontRadBarDx(idx))*(XpontRadBarDx(idx) - XpontRadDx(idx))*2);
                end
            end
        end
        %% PARTE AGGIUNTIVA PM 2 ====================================================
        %        Var = ones(1,nlay);
        %        rG_PM = zeros(1,nlay);
        %        t = 0;
        %        while abs((sum(rG_PM) - sum(Var))) > 1e-3
        %            t = t + 1;
        %            if t >= 2
        %                Var = rG_PM;
        %            end
        for i = 1 : nlay
            tmp_theta(i) = atan2(YpontRadBarDx(i),x0-XpontRadBarDx(i));
            thetaArco(i) = theta(i) - tmp_theta(i);
        end
        area_barr_withoutPM = AreaBarr;
        
        tempBr = mean(unique(mat.LayerMag.Br));  % c'è magnete o aria?? (gp .. temporaneo ott 16, 2015 )
        
        if tempBr > 0 && tempBr <= 1.25
            area_barr_withPM = area_barr_withoutPM.*mat.LayerMag.Br/1.25; % geo.Br_commercial is the Br of the real magnet
            for i = 1 : 2: length(r_all)
                idx = ceil(i/2);
                rPM = (2/3*sin(0.5*thetaArco(idx))/(0.5*thetaArco(idx))*(r_all(2*idx)^3 - r_all(2*idx - 1)^3)/(r_all(2*idx)^2 - r_all(2*idx - 1)^2));
                xG = x0 - rPM*cos(thetaArco(idx)*0.5 + tmp_theta(idx));
                if area_barr_withPM(idx) <= (r_all(2*idx)^2 - r_all(2*idx - 1)^2)*thetaArco(idx)
                    theta_tmp = area_barr_withPM(idx)/(r_all(2*idx)^2 - r_all(2*idx - 1)^2);
                    rPM = (2/3*sin(0.5*theta_tmp)/(0.5*theta_tmp)*(r_all(2*idx)^3 - r_all(2*idx - 1)^3)/(r_all(2*idx)^2 - r_all(2*idx - 1)^2));
                    rG_PM(idx) = x0 - rPM*cos(theta_tmp*0.5 + tmp_theta(idx));
                elseif area_barr_withPM(idx) <= ((r_all(2*idx)^2 - r_all(2*idx - 1)^2)*thetaArco(idx) + pi*r_semi(idx)^2)
                    x_circle = x0 - mean([X4(idx) X3(idx)]);
                    y_circle = mean([Y4(idx) Y3(idx)]);
                    theta_tmp = atan2(y_circle,x_circle);
                    theta_semiCircle = 2*(0.5*area_barr_withPM(idx) - 0.5*(r_all(2*idx)^2 - r_all(2*idx - 1)^2)*thetaArco(idx))/r_semi(idx)^2;
                    xG_semiCircleR = 4/3*r_semi(idx)*sin(theta_semiCircle/2)/theta_semiCircle;
                    rG_semiCircleA = sqrt(xG_semiCircleR^2 + (x_circle^2 + y_circle^2));
                    theta_tmp1 = asin(xG_semiCircleR/rG_semiCircleA);
                    xG_semiCircle = x0 - rG_semiCircleA*cos(theta_tmp1 + theta_tmp);
                    rG_PM(idx) = (xG*(r_all(2*idx)^2 - r_all(2*idx - 1)^2)*thetaArco(idx)*0.5 + xG_semiCircle*0.5*theta_semiCircle*r_semi(idx)^2)/(0.5*area_barr_withPM(idx));
                else
                    x_circle = x0 - mean([X4(idx) X3(idx)]);
                    y_circle = mean([Y4(idx) Y3(idx)]);
                    theta_tmp = atan2(y_circle,x_circle);
                    xG_semiCircleR = 4/3*r_semi(idx)/pi;
                    rG_semiCircleA = sqrt(xG_semiCircleR^2 + (x_circle^2 + y_circle^2));
                    theta_tmp1 = asin(xG_semiCircleR/rG_semiCircleA);
                    xG_semiCircle = x0 - rG_semiCircleA*cos(theta_tmp1 + theta_tmp);
                    AreaRes = (0.5*area_barr_withPM(idx) - 0.5*(r_all(2*idx)^2 - r_all(2*idx - 1)^2)*thetaArco(idx) - 0.5*pi*r_semi(idx)^2);
                    rG_PM(idx) = (xG*(r_all(2*idx)^2 - r_all(2*idx - 1)^2)*thetaArco(idx)*0.5 + xG_semiCircle*0.5*pi*r_semi(idx)^2 + mean([XpontRadDx XpontRadSx])*AreaRes)/(0.5*area_barr_withPM(idx));
                end
            end
            for ii = 1 : nlay
                rG(ii) = (Afe(ii)*rhoFE*rG(ii) + 0.5*area_barr_withPM(ii)*rhoPM*rG_PM(ii)) / (Afe(ii)*rhoFE + 0.5*area_barr_withPM(ii)*rhoPM);  % baricentri
            end
            M_Fe = 2*Afe*l * 1e-9 * rhoFE;   % massa ferro appeso ai ponticelli
            A_PM = cumsum(area_barr_withPM);
            M_PM = A_PM*l* 1e-9 *rhoPM;
            F_centrifuga = (M_Fe + M_PM) .* rG/1000 *  (nmax * pi/30)^2;
            
            pont = F_centrifuga/(sigma_max * l);    % mm
            for jj = 1 : 2 : length(r_all)
                idx = ceil(jj/2);
                if (pont(idx) < geo.pont0) % non disegno i ponticelli radiali il cui spessore è minore della tolleranza di lavorazione
                    pont(idx) = 0;
                end
                Counter = 0;
                while (area_barr_withoutPM(idx) < area_barr_withPM(idx) + pont(idx)*hc(idx))
                    Counter = Counter + 1;
                    Br = 0.99*mat.LayerMag.Br;
                    area_barr_withPM = area_barr_withoutPM*Br/1.25; %geo.Br_commercial is the Br of the real magnet
                    rPM = (2/3*sin(0.5*thetaArco(idx))/(0.5*thetaArco(idx))*(r_all(2*idx)^3 - r_all(2*idx - 1)^3)/(r_all(2*idx)^2 - r_all(2*idx - 1)^2));
                    xG = x0 - rPM*cos(thetaArco(idx)*0.5 + tmp_theta(idx));
                    if area_barr_withPM(idx) <= (r_all(2*idx)^2 - r_all(2*idx - 1)^2)*thetaArco(idx)
                        theta_tmp = area_barr_withPM(idx)/(r_all(2*idx)^2 - r_all(2*idx - 1)^2);
                        rPM = (2/3*sin(0.5*theta_tmp)/(0.5*theta_tmp)*(r_all(2*idx)^3 - r_all(2*idx - 1)^3)/(r_all(2*idx)^2 - r_all(2*idx - 1)^2));
                        rG_PM(idx) = x0 - rPM*cos(theta_tmp*0.5 + tmp_theta(idx));
                    elseif area_barr_withPM(idx) <= ((r_all(2*idx)^2 - r_all(2*idx - 1)^2)*thetaArco(idx) + pi*r_semi(idx)^2)
                        x_circle = x0 - mean([X4(idx) X3(idx)]);
                        y_circle = mean([Y4(idx) Y3(idx)]);
                        theta_tmp = atan2(y_circle,x_circle);
                        theta_semiCircle = 2*(0.5*area_barr_withPM(idx) - 0.5*(r_all(2*idx)^2 - r_all(2*idx - 1)^2)*thetaArco(idx))/r_semi(idx)^2;
                        xG_semiCircleR = 4/3*r_semi(idx)*sin(theta_semiCircle/2)/theta_semiCircle;
                        rG_semiCircleA = sqrt(xG_semiCircleR^2 + (x_circle^2 + y_circle^2));
                        theta_tmp1 = asin(xG_semiCircleR/rG_semiCircleA);
                        xG_semiCircle = x0 - rG_semiCircleA*cos(theta_tmp1 + theta_tmp);
                        rG_PM(idx) = (xG*(r_all(2*idx)^2 - r_all(2*idx - 1)^2)*thetaArco(idx)*0.5 + xG_semiCircle*0.5*theta_semiCircle*r_semi(idx)^2)/(0.5*area_barr_withPM(idx));
                    else
                        x_circle = x0 - mean([X4(idx) X3(idx)]);
                        y_circle = mean([Y4(idx) Y3(idx)]);
                        theta_tmp = atan2(y_circle,x_circle);
                        xG_semiCircleR = 4/3*r_semi(idx)/pi;
                        rG_semiCircleA = sqrt(xG_semiCircleR^2 + (x_circle^2 + y_circle^2));
                        theta_tmp1 = asin(xG_semiCircleR/rG_semiCircleA);
                        xG_semiCircle = x0 - rG_semiCircleA*cos(theta_tmp1 + theta_tmp);
                        AreaRes = (0.5*area_barr_withPM(idx) - 0.5*(r_all(2*idx)^2 - r_all(2*idx - 1)^2)*thetaArco(idx) - 0.5*pi*r_semi(idx)^2);
                        rG_PM(idx) = (xG*(r_all(2*idx)^2 - r_all(2*idx - 1)^2)*thetaArco(idx)*0.5 + xG_semiCircle*0.5*pi*r_semi(idx)^2 + mean([XpontRadDx XpontRadSx])*AreaRes)/(0.5*area_barr_withPM(idx));
                    end
                    if Counter == 10000
                        disp('Error');
                        break
                    end
                end
            end
            for ii = 1 : nlay
                rG(ii) = (Afe(ii)*rhoFE*rG(ii) + 0.5*area_barr_withPM(ii)*rhoPM*rG_PM(ii)) / (Afe(ii)*rhoFE + 0.5*area_barr_withPM(ii)*rhoPM);  % baricentri
            end
            M_Fe = 2*Afe*l * 1e-9 * rhoFE;   % massa ferro appeso ai ponticelli
            A_PM = cumsum(area_barr_withPM);
            M_PM = A_PM*l* 1e-9 *rhoPM;
            F_centrifuga = (M_Fe + M_PM) .* rG/1000 *  (nmax * pi/30)^2;
            if geo.radial_ribs_eval == 0
                pont = F_centrifuga/(sigma_max * l);    % mm
            else
                pont = geo.pont;
            end
            for jj = 1 : nlay
                if (pont(jj) < pont0) % non disegno i ponticelli radiali il cui spessore è minore della tolleranza di lavorazione per gli altri tipi di rotore
                    pont(jj) = 0;
                end
            end
            
            for jj = 1 : 2 : length(r_all)
                % DISEGNO DEI PONTICELLI RADIALI
                % keyboard
                x10 = x0 - r_all(jj) - racc_pont;                           % Coordinata x Nodo 10 ((Ponticello radiale_vertice in alto a dx)=(Punto
                % in basso a dx del raccordo))
                x11 = x0 - r_all(jj+1) + racc_pont;                         % Coordinata x Nodo 11 ((Ponticello radiale_vertice in alto a sx)=(Punto
                % in basso a sx del raccordo))
                if (x11 < x10) && (pont(ceil(jj/2))>0)                      % Se i due punti non si incrociano (raccordo non troppo grande rispetto alla
                    %                     % larghezza della barriera), e se lo spessore estratto dall'algoritmo per il
                    %                     % ponticello non è troppo piccolo, allora procedo al disegno del ponticello
                    %                     lpont = x10 - x11;                                      % lpont=lunghezza ponticello=larghezza barriera-2*racc_pont
                    %                     hpont = pont(ceil(jj/2));
                    %                     y10 = hpont/2;                                          % Coordinata y Nodo 10
                    %                     y11 = y10;                                              % Coordinata y Nodo 11
                    %                     XpontRadDx(ceil(jj/2)) = x10; YpontRadDx(ceil(jj/2)) = y10;
                    %                     XpontRadSx(ceil(jj/2)) = x11; YpontRadSx(ceil(jj/2)) = y10;
                    %                     x12 = x10 + 1.5 * racc_pont; y12 = y10+1.5*racc_pont;   % Punto 12 (serve solo per il disegno del raccordo: verrà cancellato)
                    %                     XpontRadBarDx(ceil(jj/2)) = XBanqdx(ceil(jj/2)); YpontRadBarDx(ceil(jj/2)) = y10 + (XpontRadBarDx(ceil(jj/2))-x10);
                    %                     x13 = x11 - 1.5 * racc_pont; y13 = y10+1.5*racc_pont;   % Punto 13 (serve solo per il disegno del raccordo: verrà cancellato)
                    %                     XpontRadBarSx(ceil(jj/2)) = XBanqsx(ceil(jj/2)); YpontRadBarSx(ceil(jj/2)) = y10 + (x11-XpontRadBarSx(ceil(jj/2)));
                    %
                    hpont = pont(ceil(jj/2));
                    y10 = hpont/2;                                          % Coordinata y Nodo 10
                    YpontRadBarDx(ceil(jj/2)) = y10 + racc_pont;             % Coordinata y pontRadBarDx
                    [x_temp,~] = intersezione_retta_circonferenza(x0,0,r_all(jj),0,YpontRadBarDx(ceil(jj/2)));
                    XpontRadBarDx(ceil(jj/2)) = 2*x0 - x_temp;
                    YpontRadBarSx(ceil(jj/2)) = y10 + racc_pont;             % Coordinata y pontRadBarDx
                    [x_temp,~] = intersezione_retta_circonferenza(x0,0,r_all(jj+1),0,YpontRadBarSx(ceil(jj/2)));
                    XpontRadBarSx(ceil(jj/2)) = 2*x0 - x_temp;
                    XpontRadDx(ceil(jj/2))  = XpontRadBarDx(ceil(jj/2))  - racc_pont;
                    YpontRadDx(ceil(jj/2))  = y10;
                    XpontRadSx(ceil(jj/2))  = XpontRadBarSx(ceil(jj/2))  + racc_pont;
                    YpontRadSx(ceil(jj/2))  = y10;
                else                                                         % Se invece i punti 10 e 11 si incrociano (ovvero raccordo troppo grande rispetto)
                    % alla larghezza barriera), non faccio il ponticello e disegno solo una linea di
                    hpont = 0;
                    XpontRadBarSx(ceil(jj/2)) = XBanqsx(ceil(jj/2));
                    YpontRadBarSx(ceil(jj/2)) = 0;
                    XpontRadBarDx(ceil(jj/2)) = XBanqdx(ceil(jj/2));
                    YpontRadBarDx(ceil(jj/2)) = 0;
                    XpontRadDx(ceil(jj/2)) = NaN;
                    YpontRadDx(ceil(jj/2)) = 0;
                    XpontRadSx(ceil(jj/2)) = NaN;
                    YpontRadSx(ceil(jj/2)) = 0;
                end
                if ceil(jj/2) == 1;
                    x = x0 - rbeta(ceil(jj/2));
                    y = hpont(1) * 0.6;
                    if y == 0
                        y = pont0;
                    end
                else
                    [xtemp,ytemp] = rot_point(-rbeta(ceil(jj/2)),0,-0.5*beta(ceil(jj/2))*pi/180);
                    x = x0 + xtemp;
                    y = ytemp;
                end
            end
        elseif mat.LayerMag.Br > 1.25
            disp('Error: wrong Br value');
        end
        
    otherwise
        
        [xrot_traf,yrot_traf]=calc_intersezione_cerchi(r, rTemp, x0);
        Dx=(r-xrot_traf(1))/5;  %per avere almeno 5 divisioni;
        xcir_plot=[r:-Dx:r*cos(pi/2/p)];
        ycir_plot=sqrt(r^2-xcir_plot.^2);
        VectCir=find(xcir_plot>=xrot_traf(1));
        x_ext_rot=xcir_plot(VectCir);
        y_ext_rot=ycir_plot(VectCir);
        
        A=[];
        for ii=1:nlay
            if ii==1
                X=[B2k(ii), XpBar2(ii), xxD2k(ii),xpont(ii),fliplr(x_ext_rot)];
                Y=[0, YpBar2(ii), yyD2k(ii),ypont(ii),fliplr(y_ext_rot)];
                %         figure(100);hold on;fill(X,Y,'r');hold off;
                A(ii)=polyarea(X,Y);
                tmp_bary(ii,:)=centroid(X',Y');         %MarcoP
                clear X Y;
            else
                X=[B1k(ii-1), XpBar1(ii-1), xxD1k(ii-1),xpont(ii-1),xrot_traf(ii-1),xrot_traf(ii),xpont(ii),xxD2k(ii),XpBar2(ii),B2k(ii)];
                Y=[0, YpBar1(ii-1), yyD1k(ii-1),ypont(ii-1),yrot_traf(ii-1),yrot_traf(ii),ypont(ii),yyD2k(ii),YpBar2(ii),0];
                %         figure(100);hold on;fill(X,Y,'r'); hold off;
                A(ii)=polyarea(X,Y);
                tmp_bary(ii,:)=centroid(X',Y');         %MarcoP
                clear X Y;
            end
        end
        Afe=cumsum(A);
        
        %MarcoP
        rG(1)=tmp_bary(1,1);     % baricentro della regione di ferro sorretta dal ponticello a I
        for ii=2:nlay
            rG(ii)=(Afe(ii-1)*rG(ii-1)+A(ii)*tmp_bary(ii,1))/Afe(ii);  % baricentri delle regioni di ferro sorrette dalle due U
        end
        
        M_Fe = 2*Afe*l * 1e-9 * rhoFE ;   % massa ferro appeso ai ponticelli
        
        F_centrifuga = M_Fe .* rG/1000 *  (nmax * pi/30)^2;
        
        if geo.radial_ribs_eval == 0
            pont = F_centrifuga/(sigma_max * l);    % mm
        else
            pont = geo.pont;
        end
        
        for jj=1:nlay
            if (pont(jj) < pont0) % non disegno i ponticelli radiali il cui spessore è minore della tolleranza di lavorazione per gli altri tipi di rotore
                pont(jj)=0;
            end
        end
        
        
        hpont=pont/2;
        rac_pont=abs(B1k-B2k)/4;
        
        for ii=1:nlay
            
            if hpont(ii)>0
                XpontRadBarSx(ii)=B1k(ii);
                YpontRadBarSx(ii)=hpont(ii)+rac_pont(ii);
                XpontRadBarDx(ii)=B2k(ii);
                YpontRadBarDx(ii)=hpont(ii)+rac_pont(ii);
                XpontRadDx(ii)=B2k(ii)-rac_pont(ii);
                YpontRadDx(ii)=hpont(ii);
                XpontRadSx(ii)=B1k(ii)+rac_pont(ii);
                YpontRadSx(ii)=hpont(ii);
                
            else
                
                
                XpontRadBarSx(ii)=B1k(ii);
                YpontRadBarSx(ii)=0;
                XpontRadBarDx(ii)=B2k(ii);
                YpontRadBarDx(ii)=0;
                XpontRadDx(ii)=NaN;
                YpontRadDx(ii)=0;
                XpontRadSx(ii)=NaN;
                YpontRadSx(ii)=0;
            end
            
        end
        
        
        %% MarcoP 07/08/2014
        % calcolo di area e baricentro delle barriere di flusso:
        for iflag=1:nlay
            if XpontRadDx(iflag)>0
                flag_PM(iflag)=1;
            else
                flag_PM(iflag)=0;
            end
        end
        % barriera ad I
        X_Ibarr=[xxD1k(1) xxD2k(1) B2k(1) B1k(1)];
        Y_Ibarr=[yyD1k(1) yyD2k(1) 0 0];
        % areaI=2*polyarea(X_Ibarr,Y_Ibarr)+pi*((B2k(1)-B1k(1))/2)^2-flag_PM(1)*(B2k(1)-B1k(1))*(XpontRadBarDx(1)-XpontRadBarSx(1));  % area totale della barriera a I
        areaI=2*polyarea(X_Ibarr,Y_Ibarr)-flag_PM(1)*(B2k(1)-B1k(1))*(XpontRadBarDx(1)-XpontRadBarSx(1));  % area della barriera a I in cui inserire i magneti (escludo i due semicerchi all'estremità della barriera)
        
        % % barriera ad U centrale
        % X_U1barr=[xxD1k(2) xxD2k(2) XpBar2(2) B2k(2) B1k(2) XpBar1(2)];
        % Y_U1barr=[yyD1k(2) yyD2k(2) YpBar2(2) 0 0 YpBar1(2)];
        % % areaU1=2*polyarea(X_U1barr,Y_U1barr)+pi*((B2k(2)-B1k(2))/2)^2-flag_PM(2)*(B2k(2)-B1k(2))*(XpontRadBarDx(2)-XpontRadBarSx(2)); % area totale della barriera a U centrale
        % areaU1=(2*YpBar2(2)*hc(2))-flag_PM(2)*(B2k(2)-B1k(2))*(XpontRadBarDx(2)-XpontRadBarSx(2)); %  area della barriera a U in cui inserire i magneti (escludo i due semicerchi all'estremità della barriera e le due "gambe" della "U")
        %
        % % barriera ad U interna
        % X_U2barr=[xxD1k(3) xxD2k(3) XpBar2(3) B2k(3) B1k(3) XpBar1(3)];
        % Y_U2barr=[yyD1k(3) yyD2k(3) YpBar2(3) 0 0 YpBar1(3)];
        % % areaU2=2*polyarea(X_U2barr,Y_U2barr)+pi*((B2k(3)-B1k(3))/2)^2-flag_PM(3)*(B2k(3)-B1k(3))*(XpontRadBarDx(3)-XpontRadBarSx(3)); % area totale della barriera a U interna
        % areaU2=(2*YpBar2(3)*hc(3))-flag_PM(3)*(B2k(3)-B1k(3))*(XpontRadBarDx(3)-XpontRadBarSx(3)); %  area della barriera a U in cui inserire i magneti (escludo i due semicerchi all'estremità della barriera e le due "gambe" della "U")
        if nlay == 1
            area_barr_withoutPM = [areaI];
        else
            % barriere ad U
            for kk=2:nlay
                areaUbars(kk-1) = (2*YpBar2(kk)*hc(kk))-flag_PM(kk)*(B2k(kk)-B1k(kk))*(XpontRadBarDx(kk)-XpontRadBarSx(kk)); %area della barriera a U in cui inserire i magneti (escludo i due semicerchi all'estremità della barriera e le due "gambe" della "U")
            end
            area_barr_withoutPM = [areaI areaUbars];
        end
        
        %         if (geo.Br>0 && geo.Br<=geo.Br_commercial)
        %             %     areaI_withPM = areaI*geo.Br/geo.Br_commercial;      %geo.Br_commercial is the Br of the real magnet
        %             %     areaU1_withPM = areaU1*geo.Br/geo.Br_commercial;    %geo.Br_commercial is the Br of the real magnet
        %             %     areaU2_withPM = areaU2*geo.Br/geo.Br_commercial;    %geo.Br_commercial is the Br of the real magnet
        %             %     area_barr_withPM = [areaI_withPM areaU1_withPM areaU2_withPM];
        %
        %             area_barr_withPM=area_barr_withoutPM*geo.Br/geo.Br_commercial;      %geo.Br_commercial is the Br of the real magnet
        %
        %             rG_PM = (B1k+B2k)/2;    % vettore con i baricentri dei magneti
        %             %ricalcolo i baricentri delle masse di ferro e magnete
        %             for ii=1:nlay
        %                 rG(ii)=(Afe(ii)*geo.rhoFE*rG(ii) + area_barr_withPM(ii)*geo.rhoPM*rG_PM(ii)) / (Afe(ii)*geo.rhoFE+area_barr_withPM(ii)*geo.rhoPM);  % baricentri delle regioni di ferro sorrette dalla I e dalle due U
        %             end
        %
        %             %% ripeto le righe 36-78 per ricalcolare i ponticelli radiali in presenza di magneti (da qui fino alla fine dello script)
        %             M_Fe = 2*Afe*l * 1e-9 * geo.rhoFE ;   % massa ferro appeso ai ponticelli
        %             A_PM=cumsum(area_barr_withPM);
        %             %     A_PM=cumsum([areaI_withPM, areaU1_withPM, areaU2_withPM]);
        %
        %             M_PM = A_PM*l* 1e-9 *geo.rhoPM;
        %             F_centrifuga = (M_Fe+M_PM) .* rG/1000 *  (nmax * pi/30)^2;
        %             sigma_max=geo.sigma_max;    %Yield strength of rotor lamination
        %
        %             pont = F_centrifuga/(sigma_max * l);    % mm
        %
        %             for jj=1:nlay
        %                 if (pont(jj) < geo.pont0) % non disegno i ponticelli radiali il cui spessore è minore della tolleranza di lavorazione
        %                     pont(jj)=0;
        %                 end
        %
        %                 % controllo per vedere se ho a disposizione area a sufficienza per
        %                 % il magnete
        %                 while (area_barr_withoutPM(jj) < area_barr_withPM(jj) + pont(jj)*hc(jj))
        %                     geo.Br = 0.99*geo.Br;
        %                     %             areaI_withPM = areaI*geo.Br/geo.Br_commercial;      %geo.Br_commercial is the Br of the real magnet
        %                     %             areaU1_withPM = areaU1*geo.Br/geo.Br_commercial;    %geo.Br_commercial is the Br of the real magnet
        %                     %             areaU2_withPM = areaU2*geo.Br/geo.Br_commercial;    %geo.Br_commercial is the Br of the real magnet
        %                     %             area_barr_withPM = [areaI_withPM areaU1_withPM areaU2_withPM];
        %                     area_barr_withPM=area_barr_withoutPM*geo.Br/geo.Br_commercial;      %geo.Br_commercial is the Br of the real magnet
        %
        %                     rG_PM = (B1k+B2k)/2;    % vettore con i baricentri dei magneti
        %                     %ricalcolo i baricentri delle masse di ferro e magnete
        %                     for ii=1:nlay
        %                         rG(ii)=(Afe(ii)*geo.rhoFE*rG(ii) + area_barr_withPM(ii)*geo.rhoPM*rG_PM(ii)) / (Afe(ii)*geo.rhoFE+area_barr_withPM(ii)*geo.rhoPM);  % baricentri delle regioni di ferro sorrette dalla I e dalle due U
        %                     end
        %
        %                     M_Fe = 2*Afe*l * 1e-9 * geo.rhoFE ;   % massa ferro appeso ai ponticelli
        %                     A_PM=cumsum(area_barr_withPM);
        %                     %             A_PM=cumsum([areaI_withPM, areaU1_withPM, areaU2_withPM]);
        %                     M_PM = A_PM*l* 1e-9 *geo.rhoPM;
        %                     F_centrifuga = (M_Fe+M_PM) .* rG/1000 *  (nmax * pi/30)^2;
        %                     sigma_max=geo.sigma_max;    %Yield strength of rotor lamination
        %
        %                     pont = F_centrifuga/(sigma_max * l);    % mm
        %                 end
        %
        %             end
        %
        %
        %             hpont=pont;
        %             rac_pont=abs(B1k-B2k)/4;
        %
        %             for ii=1:nlay
        %
        %                 if hpont(ii)>0
        %                     XpontRadBarSx(ii)=B1k(ii);
        %                     YpontRadBarSx(ii)=hpont(ii)+rac_pont(ii);
        %                     XpontRadBarDx(ii)=B2k(ii);
        %                     YpontRadBarDx(ii)=hpont(ii)+rac_pont(ii);
        %                     XpontRadDx(ii)=B2k(ii)-rac_pont(ii);
        %                     YpontRadDx(ii)=hpont(ii);
        %                     XpontRadSx(ii)=B1k(ii)+rac_pont(ii);
        %                     YpontRadSx(ii)=hpont(ii);
        %
        %                 else
        %
        %
        %                     XpontRadBarSx(ii)=B1k(ii);
        %                     YpontRadBarSx(ii)=0;
        %                     XpontRadBarDx(ii)=B2k(ii);
        %                     YpontRadBarDx(ii)=0;
        %                     XpontRadDx(ii)=NaN;
        %                     YpontRadDx(ii)=0;
        %                     XpontRadSx(ii)=NaN;
        %                     YpontRadSx(ii)=0;
        %                 end
        %
        %             end
        %             %%
        %         elseif geo.Br>geo.Br_commercial
        %             disp('Error: wrong Br value')
        %
        %         end
end

