% central point and direction of magnetization

% assign the condition for flux barrier division

% additional point
for kk=1:nlay
    if (yyD2k(kk)==YpBar2(kk))
        [a,b,c]=retta_per_2pti(XpBar2(kk),YpBar2(kk),xpont(kk),ypont(kk));
    else
        [a,b,c]=retta_per_2pti(XpBar2(kk),YpBar2(kk),xxD2k(kk),yyD2k(kk));
    end
    a2=b;
    b2=-a;
    c2=a*YpBar2(kk)-b*XpBar2(kk);
    [a3,b3,c3]=retta_per_2pti(XpBar1(kk),YpBar1(kk),xxD1k(kk),yyD1k(kk));
    [x,y]=intersezione_tra_rette(a2,b2,c2,a3,b3,c3);
    if x<XpBar1(kk)
        
        XpMag2B1(kk)=XpBar1(kk);
        YpMag2B1(kk)=YpBar1(kk);
        
    else
        XpMag2B1(kk)=x;
        YpMag2B1(kk)=y;
    end
    % internal barrier side intersection
    [a4,b4,c4]=retta_per_2pti(B2k(kk),0,XpBar2(kk),YpBar2(kk));
    [a5,b5,c5]=retta_per_2pti(B1k(kk),0,XpBar1(kk),YpBar1(kk));
    m6=b4/a4; d6=YpBar2(kk)-m6*XpBar2(kk);
    a6=-m6;
    b6=1;
    c6=-d6;
    [x,y]=intersezione_tra_rette(a5,b5,c5,a6,b6,c6);
    if (abs((YpBar1(kk)-y)+1j*(XpBar1(kk)-x))<0.5)
        XpMag1B1(kk)=XpBar1(kk);
        YpMag1B1(kk)=YpBar1(kk);
    elseif (y>YpBar1(kk))
        XpMag1B1(kk)=XpBar1(kk);
        YpMag1B1(kk)=YpBar1(kk);
    else
        %         if kk == 1
        %             XpMag1B1(1)=B2k(1);
        %             YpMag1B1(1)=YpBar1(1);
        %         else
        XpMag1B1(kk)=x;
        YpMag1B1(kk)=y;
        %         end
    end
    %     xv=70:0.1:80;
    %     yv4=-a4/b4*xv-c4/b4;
    %     yv5=-a5/b5*xv-c5/b5;
    %     yv6b=-a6/b6*xv-c6/b6;
    %     m4=-a4/b4;
    %     m6=-1/m4; d6=YpBar2(kk)-m6*XpBar2(kk);
    %     yv6=m6*xv+d6
    %     figure(100);plot(B2k(kk),0,'rs'); hold on;
    %     plot(XpBar2(kk),YpBar2(kk),'rs');
    %     plot(xv,yv4); plot(xv,yv5); plot(xv,yv6);plot(xv,yv6b,'r');
    %     plot(x,y,'cs');
    %     hold off; axis equal
end
clear a b c a1 b1 c1 a2 b2 c2 a3 b3 c3 a4 b4 c4 x y

% central point and direction of magnetization

% assign the condiction for flux barrier division

% addtitional point

xc=[];
yc=[];      % (xc,yc)=PMs centers
xair=[];
yair=[];    % (xair,yair)=rotor air zone centers
xmag=[];
ymag=[];    % (xmag,ymag)=PMs magnetization direction
Br = [];
xmagair=[];
ymagair=[]; % (xmagair,ymagair)= air zone magnetization direction (maybe useful for FEMM definition???)

fluxPortion=3;
for kk=1:nlay
    % PMs of the central segment
    xmedBar1=(B1k(kk)+B2k(kk))/2;
    ymedBar1=(temp.yvert2pt2(kk)+YpontRadBarDx(kk))/2;
    xc=[xc,xmedBar1];
    yc=[yc,ymedBar1];
    Br = [Br mat.LayerMag.Br(kk)];    % 1 block
    xmag=[xmag,cos(0)];
    ymag=[ymag,sin(0)];
    
    % Air zone at the end of the barrier (only if material~=Air and
    % Splitted internal ribs)
    if ((~strcmp(mat.LayerMag.MatName,'Air'))||(~isnan(XpontSplitDx(1,kk))))
        xtmp=(xxD1k(kk)+xxD2k(kk)+xpont(kk))/3;
        ytmp=(yyD1k(kk)+yyD2k(kk)+ypont(kk))/3;
        xair=[xair,xtmp];
        yair=[yair,ytmp];
    end
    %Br=[Br mat.LayerMag.Br(kk)]; % I'm not so sure...
    
    % Internal Radial Rib (only if material~=Air)
    if ((~strcmp(mat.LayerMag.MatName,'Air'))&&(~isnan(XpontRadDx(kk))))
        xtmp=(B1k(kk)+B2k(kk))/2;
        ytmp=(YpontRadBarSx(kk)+YpontRadSx(kk))/2;
        xair=[xair,xtmp];
        yair=[yair,ytmp];
    end
    
    % Splitted Internal Radial Rib (only if material~=Air)
    if (~strcmp(mat.LayerMag.MatName,'Air')&&~isnan(XpontSplitDx(1,kk)))
        xtmp=(B1k(kk)+B2k(kk))/2;
        ytmp=(YpontSplitBarSx(2,kk)+YpontSplitSx(2,kk))/2;
        xair=[xair,xtmp];
        yair=[yair,ytmp];
    end
    
    % Flux barrier arms (only if kk>1 && material~=Air)
    if ((kk>1)&&(~strcmp(mat.LayerMag.MatName,'Air'))&&(geo.Areavert(kk)>0))
        % PMs
        d2221 = calc_distanza_punti([XpBar2(kk),YpBar2(kk)],[xxD2k(kk),yyD2k(kk)]);
        if (d2221>0.5)
            [a_22tmp,b_22tmp,c_22tmp]=retta_per_2pti(temp.xob2pt1(kk),temp.yob2pt1(kk),temp.xob1pt2(kk),temp.yob1pt2(kk));
            [a_21tmp,b_21tmp,c_21tmp]=retta_per_2pti(temp.xob1pt1(kk),temp.yob1pt1(kk),temp.xob2pt2(kk),temp.yob2pt2(kk));
        else
            [a_22tmp,b_22tmp,c_22tmp]=retta_per_2pti(temp.xob2pt1(kk),temp.yob2pt1(kk),temp.xob1pt2(kk),temp.yob1pt2(kk));
            [a_21tmp,b_21tmp,c_21tmp]=retta_per_2pti(temp.xob1pt1(kk),temp.yob1pt1(kk),xpont(kk),ypont(kk));

        end
        [xmedBar3,ymedBar3]=intersezione_tra_rette(a_22tmp,b_22tmp,c_22tmp,a_21tmp,b_21tmp,c_21tmp);
        xc=[xc,xmedBar3];
        yc=[yc,ymedBar3];
        xmag=[xmag,cos(atan(mOrto))];
        ymag=[ymag,sin(atan(mOrto))];
        Br = [Br mat.LayerMag.Br(kk)];    % add another blocks for layers > 1
        
        % Air zone in the corner of the flux barrier
        if (YpBar2(kk)==yyD2k(kk))
            [a3,b3,c3]=retta_per_2pti(XpBar2(kk),YpBar2(kk),xpont(kk),ypont(kk));
        else
            [a3,b3,c3]=retta_per_2pti(XpBar2(kk),YpBar2(kk),xxD2k(kk),yyD2k(kk));
        end
        mOrto= b3/a3;

        xmedBar2=(XpBar1(kk)+XpBar2(kk))/2;
        ymedBar2=(YpMag1B1(kk)+YpBar1(kk))/2;
        [a2,b2,c2]=retta_per_2pti(XpMag1B1(kk),YpMag1B1(kk),XpBar2(kk),YpBar2(kk));
        mParal=-(a2/b2);
        %if strcmp (mat.LayerMag.MatName,'Bonded-Magnet') || strcmp(mat.LayerMag.MatName,'Air')
        xair=[xair,xmedBar2];
        yair=[yair,ymedBar2];
        xmagair=[xmagair,cos(atan(mParal))];
        ymagair=[ymagair,sin(atan(mParal))];
        Br = [Br mat.LayerMag.Br(kk)];
    end
end
zmag=zeros(1,size(xmag,2));

%% Cross section drawing

% figure(100);hold on;
% % plot(xo',yo','--r','LineWidth',2); axis([0 rlim 0 rlim]); axis square
% plot(B1k,0,'ob');plot(B2k,0,'ob');
% plot(xpont,ypont,'*c');
% plot(XpBar2,YpBar2,'bs');
% plot(XpBar1,YpBar1,'bs');
% plot(xTraf1,yTraf1,'*m');
% plot(xTraf2,yTraf2,'*m');
% plot(xD1k,yD1k,'ob');
% plot(XpMag2B1,YpMag2B1,'dc');
% %
% for ii=1:nlay
%    plot([XpBar1(ii),xxD1k(ii)],[YpBar1(ii),yyD1k(ii)],'b','LineWidth',2);
%    plot([XpBar2(ii),xxD2k(ii)],[YpBar2(ii),yyD2k(ii)],'b','LineWidth',2);
%    plot([B1k(ii),XpBar1(ii)],[0,YpBar1(ii)],'b','LineWidth',2);
%    plot([B2k(ii),XpBar2(ii)],[0,YpBar2(ii)],'b','LineWidth',2);
%
% end
% hold off;