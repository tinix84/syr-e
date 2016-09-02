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
    if (abs((YpBar1(kk)-y)+1j*(XpBar1(kk)-x))<0.5);
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

% evaluation of the direction of magnetization
xc=[];
yc=[];
xmag=[];
ymag=[]; Br = [];
fluxPortion=3;
for kk=1:nlay
    % #1 magnet portion per barrier
    
    [a1b,b1b,c1b]=retta_per_2pti(XpontRadBarSx(kk),YpontRadBarSx(kk),XpBar2(kk),YpBar2(kk));
    [a2b,b2b,c2b]=retta_per_2pti(XpontRadBarDx(kk),YpontRadBarDx(kk),XpMag1B1(kk),YpMag1B1(kk));
    
    [xmedBar1,ymedBar1]=intersezione_tra_rette(a1b,b1b,c1b,a2b,b2b,c2b);
    xc=[xc,xmedBar1];
    yc=[yc,ymedBar1];
    Br = [Br mat.LayerMag.Br(kk)];    % 1 block
    
    [a1,b1,c1]=retta_per_2pti(B2k(kk),0,XpBar2(kk),YpBar2(kk));
    if (DTrasl==0)
        mOrto=0;
    else
        mOrto=b1/a1;
    end
    xmag=[xmag,cos(atan(mOrto))];
    ymag=[ymag,sin(atan(mOrto))];
    %%  %%%%%%%%%%%%%%%%%%
    %% #2
    %%  %%%%%%%%%%%%%%%%%%
    %     ypMed22=(YpMag2B1(kk)+YpBar2(kk))/2;
    %     xpMed22=(XpMag2B1(kk)+XpBar2(kk))/2;
    %     %
    %     ypMed21=(YpMag1B1(kk)+YpBar2(kk))/2;
    %     xpMed21=(XpMag1B1(kk)+XpBar2(kk))/2;
    %     %
    %     [a1b,b1b,c1b]=retta_per_2pti(XpMag2B1(kk),YpMag2B1(kk),xpMed21,ypMed21);
    %     [a2b,b2b,c2b]=retta_per_2pti(XpMag1B1(kk),YpMag1B1(kk),xpMed22,ypMed22);
    %     %
    %         [xmedBar2,ymedBar2]=intersezione_tra_rette(a1b,b1b,c1b,a2b,b2b,c2b);
    %     %
    %     xc=[xc,xmedBar2];
    %     yc=[yc,ymedBar2];
    %     if (YpBar2(kk)==yyD2k(kk))
    %         [a2,b2,c2]=retta_per_2pti(XpBar2(kk),YpBar2(kk),xpont(kk),ypont(kk));
    %         mOrto=b2/a2;
    %     else
    %         [a2,b2,c2]=retta_per_2pti(XpBar2(kk),YpBar2(kk),XpMag2B1(kk),YpMag2B1(kk));
    %         mOrto=-a2/b2/2;
    %     end
    %     xmag=[xmag,cos(atan(mOrto))];
    %     ymag=[ymag,sin(atan(mOrto))];
    %%  %%%%%%%%%%%%%%%%%%%%
    %% #3
    %%  %%%%%%%%%%%%%%%%%%%%
    if (YpBar2(kk)==yyD2k(kk))
        [a3,b3,c3]=retta_per_2pti(XpBar2(kk),YpBar2(kk),xpont(kk),ypont(kk));
    else
        [a3,b3,c3]=retta_per_2pti(XpBar2(kk),YpBar2(kk),xxD2k(kk),yyD2k(kk));
    end
    mOrto=b3/a3;
    xmedBar3=(XpMag2B1(kk)+xxD2k(kk))/2;
    ymedBar3=(yyD1k(kk)+YpBar2(kk))/2;
    xc=[xc,xmedBar3];
    yc=[yc,ymedBar3];
    xmag=[xmag,cos(atan(mOrto))];
    ymag=[ymag,sin(atan(mOrto))];
    Br = [Br mat.LayerMag.Br(kk)];    % add another blocks to all layers 
    
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
% plot(XpMag1B1,YpMag1B1,'dr');
% %
% for ii=1:nlay
%    plot([XpBar1(ii),xxD1k(ii)],[YpBar1(ii),yyD1k(ii)],'b','LineWidth',2);
%    plot([XpBar2(ii),xxD2k(ii)],[YpBar2(ii),yyD2k(ii)],'b','LineWidth',2);
%    plot([B1k(ii),XpBar1(ii)],[0,YpBar1(ii)],'b','LineWidth',2);
%    plot([B2k(ii),XpBar2(ii)],[0,YpBar2(ii)],'b','LineWidth',2);
%
% end
% hold off;axis equal;
% keyboard