xc=[];
yc=[];
xair=[];
yair=[];
xmag=[];
ymag=[]; Br = [];
xmagair=[];
ymagair=[];
fluxPortion=3;
for kk=1:nlay
    % #1 magnet portion per barrier
    xmedBar1=(B1k(kk)+B2k(kk))/2;
    ymedBar1=(temp.yvert2pt2(kk)+YpontRadBarDx(kk))/2;
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
        
        
        xmedBar2=(XpBar1(kk)+XpBar2(kk))/2;
        ymedBar2=(YpMag1B1(kk)+YpBar1(kk))/2;
        [a2,b2,c2]=retta_per_2pti(XpMag1B1(kk),YpMag1B1(kk),XpBar2(kk),YpBar2(kk));
        mParal=-(a2/b2);
        %if strcmp (mat.LayerMag.MatName,'Bonded-Magnet') || strcmp(mat.LayerMag.MatName,'Air')
        if strcmp(mat.LayerMag.MatName,'Air')
            
            xc=[xc,xmedBar2];
            yc=[yc,ymedBar2];
            xmag=[xmag,cos(atan(mParal))];
            ymag=[ymag,sin(atan(mParal))];
            Br = [Br mat.LayerMag.Br(kk)];
        else
            if~(geo.Areavert(kk) ==0)
                xair=[xair,xmedBar2];
                yair=[yair,ymedBar2];
                xmagair=[xmagair,cos(atan(mParal))];
                ymagair=[ymagair,sin(atan(mParal))];
                Br = [Br mat.LayerMag.Br(kk)];
            end
        end
        
        d2221 = calc_distanza_punti([XpBar2(kk),YpBar2(kk)],[xxD2k(kk),yyD2k(kk)]);
        if (d2221>0.5)
            [a_22tmp,b_22tmp,c_22tmp]=retta_per_2pti(temp.xob2pt1(kk),temp.yob2pt1(kk),temp.xob1pt2(kk),temp.yob1pt2(kk));
            [a_21tmp,b_21tmp,c_21tmp]=retta_per_2pti(temp.xob1pt1(kk),temp.yob1pt1(kk),temp.xob2pt2(kk),temp.yob2pt2(kk));
        else
            [a_22tmp,b_22tmp,c_22tmp]=retta_per_2pti(temp.xob2pt1(kk),temp.yob2pt1(kk),temp.xob1pt2(kk),temp.yob1pt2(kk));
            [a_21tmp,b_21tmp,c_21tmp]=retta_per_2pti(temp.xob1pt1(kk),temp.yob1pt1(kk),xpont(kk),ypont(kk));
            
        end
        if ~(geo.Areavert(kk) == 0)
            [xmedBar3,ymedBar3]=intersezione_tra_rette(a_22tmp,b_22tmp,c_22tmp,a_21tmp,b_21tmp,c_21tmp);
            
            xc=[xc,xmedBar3];
            yc=[yc,ymedBar3];
            xmag=[xmag,cos(atan(mOrto))];
            ymag=[ymag,sin(atan(mOrto))];
            
            Br = [Br mat.LayerMag.Br(kk)];    % add another blocks for layers > 1
            
        end
        xmedBar4=(temp.xob1pt2(kk)+temp.xob2pt2(kk)+xpont(kk))/3;
        ymedBar4=(temp.yob1pt2(kk)+temp.yob2pt2(kk)+ypont(kk))/3;
        %if strcmp (mat.LayerMag.MatName,'Bonded-Magnet') || strcmp(mat.LayerMag.MatName,'Air')
        if strcmp(mat.LayerMag.MatName,'Air')
            xc=[xc,xmedBar4];
            yc=[yc,ymedBar4];
            xmag=[xmag,cos(0)];
            ymag=[ymag,sin(0)];
            Br = [Br mat.LayerMag.Br(kk)];
        else
            xair=[xair,xmedBar4];
            yair=[yair,ymedBar4];
            xmagair=[xmagair,cos(0)];
            ymagair=[ymagair,sin(0)];
            Br = [Br mat.LayerMag.Br(kk)];
        end
    else
        xmedBar4=(temp.xpont(kk)+temp.XpMag1B1(kk)+temp.XpBar2(kk))/3;
        ymedBar4=(temp.ypont(kk)+temp.YpMag1B1(kk)+temp.YpBar2(kk))/3;
        
        %if strcmp (mat.LayerMag.MatName,'Bonded-Magnet') || strcmp(mat.LayerMag.MatName,'Air')
        if strcmp(mat.LayerMag.MatName,'Air')
            xc=[xc,xmedBar4];
            yc=[yc,ymedBar4];
            xmag=[xmag,cos(0)];
            ymag=[ymag,sin(0)];
            Br = [Br mat.LayerMag.Br(kk)];
        else
            xair=[xair,xmedBar4];
            yair=[yair,ymedBar4];
            xmagair=[xmagair,cos(0)];
            ymagair=[ymagair,sin(0)];
            Br = [Br mat.LayerMag.Br(kk)];
        end
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
