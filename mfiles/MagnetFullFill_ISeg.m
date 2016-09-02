% central point and direction of magnetization

% assign the condiction for flux barrier division

% addtitional point

XpMag1B1=XpontRadBarSx;
YpMag1B1=YpBar2;
xc=[];
yc=[];
xmag=[];
ymag=[]; Br = [];
fluxPortion=3;
for kk=1:nlay
    % #1 magnet portion per barrier
    xmedBar1=(B1k(kk)+B2k(kk))/2;
%     if (YpontRadSx(kk)==0)
%         xc=[xc,xmedBar1];
%         yc=[yc,0];
%     else
        ymedBar1=(YpBar2(kk)+YpontRadBarDx(kk))/2;  
        xc=[xc,xmedBar1];
        yc=[yc,ymedBar1];
        Br = [Br mat.LayerMag.Br(kk)];    % 1 block 
%     end
    xmag=[xmag,cos(0)];
    ymag=[ymag,sin(0)];
    if (kk>1)
        % #2
        if (YpBar2(kk)==yyD2k(kk))
            [a3,b3,c3]=retta_per_2pti(XpBar2(kk),YpBar2(kk),xpont(kk),ypont(kk));
        else
            [a3,b3,c3]=retta_per_2pti(XpBar2(kk),YpBar2(kk),xxD2k(kk),yyD2k(kk));
        end
        mOrto=b3/a3;
        
        d2221 = calc_distanza_punti([XpBar2(kk),YpBar2(kk)],[xxD2k(kk),yyD2k(kk)]);
        if (d2221>0.5)
            [a_22tmp,b_22tmp,c_22tmp]=retta_per_2pti(XpBar2(kk),YpBar2(kk),xxD1k(kk),yyD1k(kk));
            [a_21tmp,b_21tmp,c_21tmp]=retta_per_2pti(XpBar1(kk),YpBar1(kk),xxD2k(kk),yyD2k(kk));
        else
            [a_22tmp,b_22tmp,c_22tmp]=retta_per_2pti(XpBar2(kk),YpBar2(kk),xxD1k(kk),yyD1k(kk));
            [a_21tmp,b_21tmp,c_21tmp]=retta_per_2pti(XpBar1(kk),YpBar1(kk),xpont(kk),ypont(kk));
            
        end
        [xmedBar3,ymedBar3]=intersezione_tra_rette(a_22tmp,b_22tmp,c_22tmp,a_21tmp,b_21tmp,c_21tmp);
        
        xc=[xc,xmedBar3];
        yc=[yc,ymedBar3];
        xmag=[xmag,cos(atan(mOrto))];
        ymag=[ymag,sin(atan(mOrto))];
        
        Br = [Br mat.LayerMag.Br(kk)];    % add another blocks for layers > 1
        
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
