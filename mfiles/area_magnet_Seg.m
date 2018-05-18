function [temp,geo] = area_magnet_Seg(temp,geo)

%Punti d'interesse area 1
XpBar1=temp.XpBar1;
YpBar1=temp.YpBar1;
xxD1k=temp.xxD1k;
yyD1k=temp.yyD1k;
xxD2k=temp.xxD2k;
yyD2k=temp.yyD2k;
XpBar2=temp.XpBar2;
YpBar2=temp.YpBar2;

%Punti d'interesse area 2
B1k=temp.B1k;
B2k=temp.B2k;
XpMag1B1=temp.XpMag1B1;
YpMag1B1=temp.YpMag1B1;

YpontRadSx    = temp.YpontRadSx;
XpontRadSx    = temp.XpontRadSx;
YpontRadDx    = temp.YpontRadDx;
XpontRadDx    = temp.XpontRadDx;
XpontRadBarDx = temp.XpontRadBarDx;
YpontRadBarDx = temp.YpontRadBarDx;
XpontRadBarSx = temp.XpontRadBarSx;
YpontRadBarSx = temp.YpontRadBarSx;

XpontSplitBarSx = temp.XpontSplitBarSx;
YpontSplitBarSx = temp.YpontSplitBarSx;
XpontSplitBarDx = temp.XpontSplitBarDx;
YpontSplitBarDx = temp.YpontSplitBarDx;
XpontSplitDx    = temp.XpontSplitDx;
YpontSplitDx    = temp.YpontSplitDx;
XpontSplitSx    = temp.XpontSplitSx;
YpontSplitSx    = temp.YpontSplitSx;

p=geo.p;
nlay=geo.nlay;
pont0=geo.pont0;
x0=geo.x0;
dob=geo.dob;
dvert=geo.dvert;

for ii= 1:nlay
    if ii==1
        
        a1=0;
        b1=0;
        c1=0;
        
    else
        
        [a1(ii),b1(ii),c1(ii)]=retta_per_2pti(XpBar1(ii),YpBar1(ii),xxD1k(ii),yyD1k(ii));
        
        m1(ii) = -(a1(ii)/b1(ii));
        q1(ii) = -(c1(ii)/b1(ii));
        m1o(ii)= -(m1(ii)).^-1;
        q1o(ii)=  (YpBar2(ii)-m1o(ii)*XpBar2(ii));
        
        xort(ii) = -((q1(ii)-q1o(ii))/(m1(ii)-m1o(ii)));
        yort(ii) =  (m1(ii)*xort(ii))+q1(ii);
        
    end
    
end


tempmatvert=[];
tempmatob=[];
Mag=[];

% I nomi delle coordinate dei punti determinati sono da intendersi riferite ai punti
% prima della rotazione.I nomi dei modificatori invece sono riferiti ai
% punti già ruotati.

for ii= 1:nlay
    if ii==1
        
        %% Determinazione punti per il tracciamento dei rettangoli
        
        % calcolo area totale per FBS
        if geo.delta_FBS~=0
            xvert1pt1(ii)=XpontRadBarSx(ii);
            xvert1pt2(ii)=XpontRadBarSx(ii);
            yvert1pt1(ii)=YpontRadBarSx(ii);
            yvert1pt2(ii)=1*(YpontSplitBarSx(2,ii)-pont0/4-YpontRadBarSx(ii))+YpontRadBarSx(ii);
            xvert2pt1(ii)=XpontRadBarDx(ii);
            xvert2pt2(ii)=XpontRadBarDx(ii);
            yvert2pt1(ii)=YpontRadBarDx(ii);
            yvert2pt2(ii)=1*(YpontSplitBarDx(2,ii)-pont0/4-YpontRadBarDx(ii))+YpontRadBarDx(ii);
            X=[xvert1pt1(ii),xvert1pt2(ii),xvert2pt2(ii),xvert2pt1(ii),xvert1pt1(ii)];
            Y=[yvert1pt1(ii),yvert1pt2(ii),yvert2pt2(ii),yvert2pt1(ii),yvert1pt1(ii)];
            AreaObMax(ii)=polyarea(X,Y);
            dob(ii)=dob(ii)*geo.Areaob0(ii)/AreaObMax(ii);
        end
        
        xvert1pt1(ii)=XpontRadBarSx(ii);
        xvert1pt2(ii)=XpontRadBarSx(ii);
        yvert1pt1(ii)=YpontRadBarSx(ii);
        yvert1pt2(ii)=dob(ii)*(YpontSplitBarSx(2,ii)-pont0/4-YpontRadBarSx(ii))+YpontRadBarSx(ii);
        xvert2pt1(ii)=XpontRadBarDx(ii);
        xvert2pt2(ii)=XpontRadBarDx(ii);
        yvert2pt1(ii)=YpontRadBarDx(ii);
        yvert2pt2(ii)=dob(ii)*(YpontSplitBarDx(2,ii)-pont0/4-YpontRadBarDx(ii))+YpontRadBarDx(ii);
        
        %calcolo aree
        %Areaob(ii,1)=polyarea([xvert1pt1(ii),xvert1pt2(ii),xvert2pt2(ii),xvert2pt1(ii),xvert1pt1(ii)],[yvert1pt1(ii),yvert1pt2(ii),yvert2pt2(ii),yvert2pt1(ii),yvert1pt1(ii)]);
        Areavert(ii,1) = 0;
        %Areatot(ii,1)=polyarea([B1k(ii),XpBar1(ii),xxD1k(ii),xxD2k(ii),XpBar2(ii),B2k(ii),B1k(ii)],[0,YpBar1(ii),yyD1k(ii),yyD2k(ii),YpBar2(ii),0,0]);
        
        Areaob(ii,1)=polyarea([xvert1pt1(ii),xvert1pt2(ii),xvert2pt2(ii),xvert2pt1(ii),xvert1pt1(ii)],[yvert1pt1(ii),yvert1pt2(ii),yvert2pt2(ii),yvert2pt1(ii),yvert1pt1(ii)]);

        if  isequal(geo.Areaob0(ii),0)
            geo.Areavert0(ii)=Areavert(ii,1);
            geo.Areaob0(ii)=Areaob(ii,1);
            geo.Areatot(ii)=[Areavert(ii,1)+Areaob(ii,1)];
        end
        
        if geo.delta_FBS==0
            dob(ii)=dob(ii)*geo.Areaob0(ii)/Areaob(ii);
            if Areavert(ii)~=0
                dvert(ii)=dvert(ii)*geo.Areavert0(ii)/Areavert(ii);
            end
        end

        geo.Areavert(ii)=Areavert(ii,1);
        geo.Areaob(ii)=Areaob(ii,1);

        %costruzione matrice dei lati
        Mag=[Mag;
            xvert1pt1(ii) yvert1pt1(ii) xvert1pt2(ii) yvert1pt2(ii) NaN NaN eps;
            xvert1pt2(ii) yvert1pt2(ii) xvert2pt2(ii) yvert2pt2(ii) NaN NaN eps;
            xvert2pt2(ii) yvert2pt2(ii) xvert2pt1(ii) yvert2pt1(ii) NaN NaN eps;
            xvert2pt1(ii) yvert2pt1(ii) xvert1pt1(ii) yvert1pt1(ii) NaN NaN eps;
            NaN NaN NaN NaN NaN NaN  eps;
            NaN NaN NaN NaN NaN NaN  eps;
            NaN NaN NaN NaN NaN NaN  eps;
            NaN NaN NaN NaN NaN NaN  eps];

        %eps completa la matrice Mag in modo tale che abbia stesso numero di colonne di rotor

        temp.xvert1pt2(ii)=xvert1pt2(ii);
        temp.yvert1pt2(ii)=yvert1pt2(ii);
        temp.xvert2pt2(ii)=xvert2pt2(ii);
        temp.yvert2pt2(ii)=yvert2pt2(ii);
        
        
        % compatibilità: aggiunta dei punti obliqui della prima barriera,
        % con valore zero
        
        temp.xob1pt1(ii)=0;
        temp.yob1pt1(ii)=0;
        temp.xob1pt2(ii)=0;
        temp.yob1pt2(ii)=0;
        temp.xob2pt1(ii)=0;
        temp.yob2pt1(ii)=0;
        temp.xob2pt2(ii)=0;
        temp.yob2pt2(ii)=0;
        
        beta1(ii)=0;
        beta2(ii)=0;
    else
        
        % calcolo area totale per FBS
        if geo.delta_FBS~=0
            % tratto centrale
            xvert1pt1(ii)=XpontRadBarSx(ii);
            xvert1pt2(ii)=XpontRadBarSx(ii);
            yvert1pt1(ii)=YpontRadBarSx(ii);
            yvert1pt2(ii)=1*(YpontSplitBarSx(2,ii)-pont0/4-YpontRadBarSx(ii))+YpontRadBarSx(ii);
            xvert2pt1(ii)=XpontRadBarDx(ii);
            xvert2pt2(ii)=XpontRadBarDx(ii);
            yvert2pt1(ii)=YpontRadBarDx(ii);
            yvert2pt2(ii)=1*(YpontSplitBarDx(2,ii)-pont0/4-YpontRadBarDx(ii))+YpontRadBarDx(ii);
            X=[xvert1pt1(ii),xvert1pt2(ii),xvert2pt2(ii),xvert2pt1(ii),xvert1pt1(ii)];
            Y=[yvert1pt1(ii),yvert1pt2(ii),yvert2pt2(ii),yvert2pt1(ii),yvert1pt1(ii)];
            AreaObMax(ii)=polyarea(X,Y);
            dob(ii)=dob(ii)*geo.Areaob0(ii)/AreaObMax(ii);
            
            % tratto esterno
            xrif(ii) = xxD2k(ii)+(0-yyD2k(ii))/m1o(ii);
            beta1(ii)=atan(yyD1k(ii)/(xrif(ii)-xxD1k(ii)));
            beta2(ii)=atan(yyD2k(ii)/(xrif(ii)-xxD2k(ii)));
            if beta1(ii)>beta2(ii)
                [a2(ii),b2(ii),c2(ii)]=retta_per_2pti(XpBar2(ii),YpBar2(ii),xxD2k(ii),yyD2k(ii));
                m2(ii) =  -(a2(ii)/b2(ii));
                q2(ii) =  -(c2(ii)/b2(ii));
                mp(ii) =  m1o(ii);
                qp(ii) =  (yyD2k(ii)-mp(ii)*xxD2k(ii));
                [a1(ii),b1(ii),c1(ii)]=retta_per_2pti(XpBar1(ii),YpBar1(ii),xxD1k(ii),yyD1k(ii));
                m1(ii) = -(a1(ii)/b1(ii));
                q1(ii) = -(c1(ii)/b1(ii));
                xp(ii) = -((q1(ii)-qp(ii))/(m1(ii)-mp(ii)));
                yp(ii) =  (m1(ii)*xp(ii))+q1(ii);
                xob1pt1(ii)=xort(ii);
                xob1pt2(ii)=xort(ii)+1*(xp(ii)-xort(ii)-pont0/4);
                yob1pt1(ii)=m1(ii)*xob1pt1(ii)+q1(ii);
                yob1pt2(ii)=m1(ii)*xob1pt2(ii)+q1(ii);
                xob2pt1(ii)=XpBar2(ii);
                xob2pt2(ii)=XpBar2(ii)+1*(xxD2k(ii)-XpBar2(ii)-pont0/4);
                yob2pt1(ii)=m2(ii)*xob2pt1(ii)+q2(ii);
                yob2pt2(ii)=m2(ii)*xob2pt2(ii)+q2(ii);
            else
                mp(ii)=  m1o(ii);
                qp(ii)= (yyD1k(ii)-mp(ii)*xxD1k(ii));
                m2(ii) =  m1(ii);
                q2(ii) = (yyD2k(ii)-m2(ii)*xxD2k(ii));
                xp(ii) = -((q2(ii)-qp(ii))/(m2(ii)-mp(ii)));
                yp(ii) =  (m2(ii)*xp(ii))+q2(ii);
                xob1pt1(ii)=xort(ii);
                xob1pt2(ii)=xort(ii)+dvert(ii)*(xxD1k(ii)-xort(ii)-pont0/4);
                yob1pt1(ii)=m1(ii)*xob1pt1(ii)+q1(ii);
                yob1pt2(ii)=m1(ii)*xob1pt2(ii)*q1(ii);
                xob2pt1(ii)=XpBar2(ii);
                xob2pt2(ii)=XpBar2(ii)+dvert(ii)*(xp(ii)-XpBar2(ii)-pont0/4);
                yob2pt1(ii)=m2(ii)*xob2pt1(ii)+q2(ii);
                yob2pt2(ii)=m2(ii)*xob2pt2(ii)+q2(ii);
            end
            X=[xob1pt1(ii),xob1pt2(ii),xob2pt2(ii),xob2pt1(ii),xob1pt1(ii)];
            Y=[yob1pt1(ii),yob1pt2(ii),yob2pt2(ii),yob2pt1(ii),yob1pt1(ii)];
            AreaVertMax(ii)=polyarea(X,Y);
            dvert(ii)=dvert(ii)*geo.Areavert0(ii)/AreaVertMax(ii);
        end
        
%         xvert1=B1k(ii);
%         yvert1=linspace(0,dob(ii)*(YpMag1B1(ii)-pont0/4-YpontRadBarSx(ii))+YpontRadBarSx(ii),k);
%         
%         xvert2=B2k(ii);
%         yvert2=linspace(0,dob(ii)*(YpBar2(ii)-pont0/4-YpontRadBarDx(ii))+YpontRadBarDx(ii),k);
%         
%         tempmatvert=[tempmatvert;
%             xvert1,yvert1,xvert2,yvert2];
%         
%         xvert1pt1(ii)=tempmatvert(ii,1);
%         xvert1pt2(ii)=tempmatvert(ii,1);
%         yvert1pt1(ii)=tempmatvert(ii,2);
%         yvert1pt2(ii)=tempmatvert(ii,k+1);
%         xvert2pt1(ii)=tempmatvert(ii,k+2);
%         xvert2pt2(ii)=tempmatvert(ii,k+2);
%         yvert2pt1(ii)=tempmatvert(ii,k+3);
%         yvert2pt2(ii)=tempmatvert(ii,end);
        
        xvert1pt1(ii)=XpontRadBarSx(ii);
        xvert1pt2(ii)=XpontRadBarSx(ii);
        yvert1pt1(ii)=YpontRadBarSx(ii);
        yvert1pt2(ii)=dob(ii)*(YpontSplitBarSx(2,ii)-pont0/4-YpontRadBarSx(ii))+YpontRadBarSx(ii);
        xvert2pt1(ii)=XpontRadBarDx(ii);
        xvert2pt2(ii)=XpontRadBarDx(ii);
        yvert2pt1(ii)=YpontRadBarDx(ii);
        yvert2pt2(ii)=dob(ii)*(YpontSplitBarDx(2,ii)-pont0/4-YpontRadBarDx(ii))+YpontRadBarDx(ii);
        
        xrif(ii) = xxD2k(ii)+(0-yyD2k(ii))/m1o(ii);
        beta1(ii)=atan(temp.yyD1k(ii)/(xrif(ii)-temp.xxD1k(ii)));
        beta2(ii)=atan(temp.yyD2k(ii)/(xrif(ii)-temp.xxD2k(ii)));
        
        if beta1(ii)>beta2(ii)
            [a2(ii),b2(ii),c2(ii)]=retta_per_2pti(XpBar2(ii),YpBar2(ii),xxD2k(ii),yyD2k(ii));
            
            m2(ii) =  -(a2(ii)/b2(ii));
            q2(ii) =  -(c2(ii)/b2(ii));
            mp(ii) =  m1o(ii);
            qp(ii) =  (yyD2k(ii)-mp(ii)*xxD2k(ii));
            
            [a1(ii),b1(ii),c1(ii)]=retta_per_2pti(XpBar1(ii),YpBar1(ii),xxD1k(ii),yyD1k(ii));
            
            m1(ii) = -(a1(ii)/b1(ii));
            q1(ii) = -(c1(ii)/b1(ii));
            
            xp(ii) = -((q1(ii)-qp(ii))/(m1(ii)-mp(ii)));
            yp(ii) =  (m1(ii)*xp(ii))+q1(ii);
            
%             xrett1 = linspace(xort(ii),(xort(ii)+dvert(ii)*((xp(ii)-xort(ii))-(pont0/4))),k);
%             yrett1 = m1(ii)*xrett1+q1(ii);
%             
%             xrett2 = linspace(XpBar2(ii),(XpBar2(ii)+dvert(ii)*(xxD2k(ii)-XpBar2(ii)-(pont0/4))),k);
%             yrett2 = m2(ii)*xrett2+q2(ii);
%             
%             tempmatob=[tempmatob;
%                 xrett1,yrett1,xrett2,yrett2];
            
            %xob1pt1(ii)=tempmatob(ii-1,1);
            %xob1pt2(ii)=tempmatob(ii-1,k);
            %yob1pt1(ii)=tempmatob(ii-1,k+1);
            %yob1pt2(ii)=tempmatob(ii-1,2*k);
            %xob2pt1(ii)=tempmatob(ii-1,2*k+1);
            %xob2pt2(ii)=tempmatob(ii-1,3*k);
            %yob2pt1(ii)=tempmatob(ii-1,3*k+1);
            %yob2pt2(ii)=tempmatob(ii-1,end);
            
            xob1pt1(ii)=xort(ii);
            xob1pt2(ii)=xort(ii)+dvert(ii)*(xp(ii)-xort(ii)-pont0/4);
            yob1pt1(ii)=m1(ii)*xob1pt1(ii)+q1(ii);
            yob1pt2(ii)=m1(ii)*xob1pt2(ii)+q1(ii);
            xob2pt1(ii)=XpBar2(ii);
            xob2pt2(ii)=XpBar2(ii)+dvert(ii)*(xxD2k(ii)-XpBar2(ii)-pont0/4);
            yob2pt1(ii)=m2(ii)*xob2pt1(ii)+q2(ii);
            yob2pt2(ii)=m2(ii)*xob2pt2(ii)+q2(ii);
            
        else
            
            mp(ii)=  m1o(ii);
            qp(ii)= (yyD1k(ii)-mp(ii)*xxD1k(ii));
            m2(ii) =  m1(ii);
            q2(ii) = (yyD2k(ii)-m2(ii)*xxD2k(ii));
            
            xp(ii) = -((q2(ii)-qp(ii))/(m2(ii)-mp(ii)));
            yp(ii) =  (m2(ii)*xp(ii))+q2(ii);
            
%             xrett1 = linspace(xort(ii),(xort(ii)+dvert(ii)*((xxD1k(ii)-xort(ii))-(pont0/4))),k);
%             yrett1 = m1(ii)*xrett1+q1(ii);
%             
%             xrett2 = linspace(XpBar2(ii),(XpBar2(ii)+dvert(ii)*((xp(ii)-XpBar2(ii))-(pont0/4))),k);
%             yrett2 = m2(ii)*xrett2+q2(ii);
%             
%             tempmatob=[tempmatob;
%                 xrett1,yrett1,xrett2,yrett2];
            
            %xob1pt1(ii)=tempmatob(ii-1,1);
            %xob1pt2(ii)=tempmatob(ii-1,k);
            %yob1pt1(ii)=tempmatob(ii-1,k+1);
            %yob1pt2(ii)=tempmatob(ii-1,2*k);
            %xob2pt1(ii)=tempmatob(ii-1,2*k+1);
            %xob2pt2(ii)=tempmatob(ii-1,3*k);
            %yob2pt1(ii)=tempmatob(ii-1,3*k+1);
            %yob2pt2(ii)=tempmatob(ii-1,end);
            
            xob1pt1(ii)=xort(ii);
            xob1pt2(ii)=xort(ii)+dvert(ii)*(xxD1k(ii)-xort(ii)-pont0/4);
            yob1pt1(ii)=m1(ii)*xob1pt1(ii)+q1(ii);
            yob1pt2(ii)=m1(ii)*xob1pt2(ii)*q1(ii);
            xob2pt1(ii)=XpBar2(ii);
            xob2pt2(ii)=XpBar2(ii)+dvert(ii)*(xp(ii)-XpBar2(ii)-pont0/4);
            yob2pt1(ii)=m2(ii)*xob2pt1(ii)+q2(ii);
            yob2pt2(ii)=m2(ii)*xob2pt2(ii)+q2(ii);
        end
        
        %% Controlli di sicurezza per compenetrazione aree
        
        if yob2pt2(ii)<YpBar2(ii)
            xob1pt1(ii)=xort(ii);
            yob1pt1(ii)=yort(ii);
            xob1pt2(ii)=xort(ii);
            yob1pt2(ii)=yort(ii);
            xob2pt1(ii)=XpBar2(ii);
            yob2pt1(ii)=YpBar2(ii);
            xob2pt2(ii)=XpBar2(ii);
            yob2pt2(ii)=YpBar2(ii);
        end
        
        if  yyD1k(ii)<yort(ii)
            xxD1k(kk)=xort(kk);
            yyD1k(kk)=yort(kk);
            xob1pt1(ii)=xort(ii);
            yob1pt1(ii)=yort(ii);
            xob1pt2(ii)=xort(ii);
            yob1pt2(ii)=yort(ii);
            xob2pt1(ii)=XpBar2(ii);
            yob2pt1(ii)=YpBar2(ii);
            xob2pt2(ii)=XpBar2(ii);
            yob2pt2(ii)=YpBar2(ii);
            
        end
        
        %% Completamento matrice Mag dei collegamenti
        
            
        Mag=[Mag;
            xvert1pt1(ii) yvert1pt1(ii) xvert1pt2(ii) yvert1pt2(ii) NaN NaN eps;
            xvert1pt2(ii) yvert1pt2(ii) xvert2pt2(ii) yvert2pt2(ii) NaN NaN eps;
            xvert2pt2(ii) yvert2pt2(ii) xvert2pt1(ii) yvert2pt1(ii) NaN NaN eps;
            xvert2pt1(ii) yvert2pt1(ii) xvert1pt1(ii) yvert1pt1(ii) NaN NaN eps];

        Areaob(ii,1)=polyarea([xvert1pt1(ii),xvert1pt2(ii),xvert2pt2(ii),xvert2pt1(ii),xvert1pt1(ii)],[yvert1pt1(ii),yvert1pt2(ii),yvert2pt2(ii),yvert2pt1(ii),yvert1pt1(ii)]);

        if dvert(ii)~=0
            Mag=[Mag;
                xob1pt1(ii) yob1pt1(ii) xob1pt2(ii) yob1pt2(ii) NaN NaN eps;
                xob1pt2(ii) yob1pt2(ii) xob2pt2(ii) yob2pt2(ii) NaN NaN eps;
                xob2pt2(ii) yob2pt2(ii) xob2pt1(ii) yob2pt1(ii) NaN NaN eps;
                xob2pt1(ii) yob2pt1(ii) xob1pt1(ii) yob1pt1(ii) NaN NaN eps];

            Areavert(ii,1)=polyarea([xob1pt1(ii),xob1pt2(ii),xob2pt2(ii),xob2pt1(ii),xob1pt1(ii)],[yob1pt1(ii),yob1pt2(ii),yob2pt2(ii),yob2pt1(ii),yob1pt1(ii)]);
        else
            Areavert(ii,1)=0;

        end

        %Areatot(ii,1)=polyarea([B1k(ii),XpBar1(ii),xxD1k(ii),xxD2k(ii),XpBar2(ii),B2k(ii),B1k(ii)],[0,YpBar1(ii),yyD1k(ii),yyD2k(ii),YpBar2(ii),0,0]);

        if  isequal(geo.Areaob0(ii),0)
            geo.Areavert0(ii)=Areavert(ii,1);
            geo.Areaob0(ii)=Areaob(ii,1);
            geo.Areatot(ii)=[Areavert(ii,1)+Areaob(ii,1)];

        end

        geo.Areavert(ii)=Areavert(ii,1);
        geo.Areaob(ii)=Areaob(ii,1);
        
        %salvataggio elementi utili in seguito
        temp.xob1pt1(ii)=xob1pt1(ii);
        temp.yob1pt1(ii)=yob1pt1(ii);
        temp.xob1pt2(ii)=xob1pt2(ii);
        temp.yob1pt2(ii)=yob1pt2(ii);
        temp.xob2pt1(ii)=xob2pt1(ii);
        temp.yob2pt1(ii)=yob2pt1(ii);
        temp.xob2pt2(ii)=xob2pt2(ii);
        temp.yob2pt2(ii)=yob2pt2(ii);
        temp.xvert1pt2(ii)=xvert1pt2(ii);
        temp.yvert1pt2(ii)=yvert1pt2(ii);
        temp.xvert2pt2(ii)=xvert2pt2(ii);
        temp.yvert2pt2(ii)=yvert2pt2(ii);
        
    end
    
    
end

% salvataggio elementi utili in seguito
temp.Mag=Mag;
temp.beta1=beta1;
temp.beta2=beta2;

end