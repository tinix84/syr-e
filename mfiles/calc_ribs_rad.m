%% 

switch geo.RotType

    case 'Circular'
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
        
        M_Fe = 2*Afe*l * 1e-9 * geo.rhoFE ;   % massa ferro appeso ai ponticelli
        
        F_centrifuga = M_Fe .* rG/1000 *  (nmax * pi/30)^2;
        
        pont = F_centrifuga/(sigma_max * l);    % mm
        
        for jj = 1 : nlay
            if (pont(jj) < pont0) % non disegno i ponticelli radiali il cui spessore è minore della tolleranza di lavorazione per gli altri tipi di rotore
                pont(jj) = 0;
            end
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
        
        M_Fe = 2*Afe*l * 1e-9 * geo.rhoFE ;   % massa ferro appeso ai ponticelli
        
        F_centrifuga = M_Fe .* rG/1000 *  (nmax * pi/30)^2;
        
        pont = F_centrifuga/(sigma_max * l);    % mm
        
        for jj=1:nlay
            if (pont(jj) < pont0) % non disegno i ponticelli radiali il cui spessore è minore della tolleranza di lavorazione per gli altri tipi di rotore
                pont(jj)=0;
            end
        end
        
        
        hpont=pont;
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
        
        % barriere ad U
        for kk=2:nlay
            areaUbars(kk-1) = (2*YpBar2(kk)*hc(kk))-flag_PM(kk)*(B2k(kk)-B1k(kk))*(XpontRadBarDx(kk)-XpontRadBarSx(kk)); %area della barriera a U in cui inserire i magneti (escludo i due semicerchi all'estremità della barriera e le due "gambe" della "U")
        end
        
        area_barr_withoutPM = [areaI areaUbars];
        
        if (geo.Br>0 && geo.Br<=geo.Br_commercial)
            %     areaI_withPM = areaI*geo.Br/geo.Br_commercial;      %geo.Br_commercial is the Br of the real magnet
            %     areaU1_withPM = areaU1*geo.Br/geo.Br_commercial;    %geo.Br_commercial is the Br of the real magnet
            %     areaU2_withPM = areaU2*geo.Br/geo.Br_commercial;    %geo.Br_commercial is the Br of the real magnet
            %     area_barr_withPM = [areaI_withPM areaU1_withPM areaU2_withPM];
            
            area_barr_withPM=area_barr_withoutPM*geo.Br/geo.Br_commercial;      %geo.Br_commercial is the Br of the real magnet
            
            rG_PM = (B1k+B2k)/2;    % vettore con i baricentri dei magneti
            %ricalcolo i baricentri delle masse di ferro e magnete
            for ii=1:nlay
                rG(ii)=(Afe(ii)*geo.rhoFE*rG(ii) + area_barr_withPM(ii)*geo.rhoPM*rG_PM(ii)) / (Afe(ii)*geo.rhoFE+area_barr_withPM(ii)*geo.rhoPM);  % baricentri delle regioni di ferro sorrette dalla I e dalle due U
            end
            
            %% ripeto le righe 36-78 per ricalcolare i ponticelli radiali in presenza di magneti (da qui fino alla fine dello script)
            M_Fe = 2*Afe*l * 1e-9 * geo.rhoFE ;   % massa ferro appeso ai ponticelli
            A_PM=cumsum(area_barr_withPM);
            %     A_PM=cumsum([areaI_withPM, areaU1_withPM, areaU2_withPM]);
            
            M_PM = A_PM*l* 1e-9 *geo.rhoPM;
            F_centrifuga = (M_Fe+M_PM) .* rG/1000 *  (nmax * pi/30)^2;
            sigma_max=geo.sigma_max;    %Yield strength of rotor lamination
            
            pont = F_centrifuga/(sigma_max * l);    % mm
            
            for jj=1:nlay
                if (pont(jj) < geo.pont0) % non disegno i ponticelli radiali il cui spessore è minore della tolleranza di lavorazione
                    pont(jj)=0;
                end
                
                % controllo per vedere se ho a disposizione area a sufficienza per
                % il magnete
                while (area_barr_withoutPM(jj) < area_barr_withPM(jj) + pont(jj)*hc(jj))
                    geo.Br = 0.99*geo.Br;
                    %             areaI_withPM = areaI*geo.Br/geo.Br_commercial;      %geo.Br_commercial is the Br of the real magnet
                    %             areaU1_withPM = areaU1*geo.Br/geo.Br_commercial;    %geo.Br_commercial is the Br of the real magnet
                    %             areaU2_withPM = areaU2*geo.Br/geo.Br_commercial;    %geo.Br_commercial is the Br of the real magnet
                    %             area_barr_withPM = [areaI_withPM areaU1_withPM areaU2_withPM];
                    area_barr_withPM=area_barr_withoutPM*geo.Br/geo.Br_commercial;      %geo.Br_commercial is the Br of the real magnet
                    
                    rG_PM = (B1k+B2k)/2;    % vettore con i baricentri dei magneti
                    %ricalcolo i baricentri delle masse di ferro e magnete
                    for ii=1:nlay
                        rG(ii)=(Afe(ii)*geo.rhoFE*rG(ii) + area_barr_withPM(ii)*geo.rhoPM*rG_PM(ii)) / (Afe(ii)*geo.rhoFE+area_barr_withPM(ii)*geo.rhoPM);  % baricentri delle regioni di ferro sorrette dalla I e dalle due U
                    end
                    
                    M_Fe = 2*Afe*l * 1e-9 * geo.rhoFE ;   % massa ferro appeso ai ponticelli
                    A_PM=cumsum(area_barr_withPM);
                    %             A_PM=cumsum([areaI_withPM, areaU1_withPM, areaU2_withPM]);
                    M_PM = A_PM*l* 1e-9 *geo.rhoPM;
                    F_centrifuga = (M_Fe+M_PM) .* rG/1000 *  (nmax * pi/30)^2;
                    sigma_max=geo.sigma_max;    %Yield strength of rotor lamination
                    
                    pont = F_centrifuga/(sigma_max * l);    % mm
                end
                
            end
            
            
            hpont=pont;
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
            %%
        elseif geo.Br>geo.Br_commercial
            disp('Error: wrong Br value')
            
        end
end

