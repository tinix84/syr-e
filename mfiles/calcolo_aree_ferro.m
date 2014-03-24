%% Calcolo delle aree di ferro di barriera

Afe=zeros(1,nlay);
xBarrier1=[];
yBarrier1=[];
xBarrier2=[];
yBarrier2=[];


for i=1:nlay
       Dangle=0.01;

    if (i==1)
        [xc2,yc2,r2,angleA2,angleB2]=circonferenza_per_3_pti(B2k(i),0,xTraf02(i),yTraf02(i),xxB2k_mean(i),yyB2k_mean(i));
        if (angleA2<angleB2)
            angleB2=-2*pi+angleB2;
        end
        xBarrier2=xc2+r2*cos([angleB2:Dangle:angleA2]);
        yBarrier2=yc2+r2*sin([angleB2:Dangle:angleA2]);
        
        for ik=1:length(xBarrier2)
            if ik==1
                Dx2(ik)=xBarrier2(1);
            else
                Dx2(ik)=xBarrier2(ik-1)-xBarrier2(ik);
            end
        end
        dim=length(Dx2); meanDx2=mean(Dx2(ceil(dim/6):dim));
        A2=sum(yBarrier2.*meanDx2);
%         A22=trapz(xBarrier2,yBarrier2);
%         0.5*(xTraf02(i)-B2k(i))*yTraf02(i)
        Afe(i)=0.5*(xr^2*atan2(yTraf02(i),xTraf02(i))-yTraf02(i)*xTraf02(i))+A2;     % Area mezza barriera di flux
        %%
    else
        [xc1,yc1,r1,angleA1,angleB1]=circonferenza_per_3_pti(B1k(i-1),0,xTraf01(i-1),yTraf01(i-1),xxB1k_mean(i-1),yyB1k_mean(i-1));
        [xc2,yc2,r2,angleA2,angleB2]=circonferenza_per_3_pti(B2k(i),0,xTraf02(i),yTraf02(i),xxB2k_mean(i),yyB2k_mean(i));
        
        if (angleA1<angleB1)
            angleB1=-2*pi+angleB1;
        end
        xBarrier1=xc1+r1*cos([angleB1:Dangle:angleA1]);
        yBarrier1=yc1+r1*sin([angleB1:Dangle:angleA1]);
        
        for ik=1:length(xBarrier1)
            if ik==1
                Dx1(ik)=xBarrier1(1);
            else
                Dx1(ik)=xBarrier1(ik)-xBarrier1(ik-1);
            end
        end
        dim=length(Dx1); meanDx1=mean(Dx1(ceil(dim/6):dim));        
        A1=sum(yBarrier1.*meanDx1);
        
        if (angleA2<angleB2)
            angleB2=-2*pi+angleB2;
        end
        xBarrier2=xc2+r2*cos([angleB2:Dangle:angleA2]);
        yBarrier2=yc2+r2*sin([angleB2:Dangle:angleA2]);
        
        for ik=1:length(xBarrier2)
            if ik==1
                Dx2(ik)=xBarrier2(1);
            else
                Dx2(ik)=xBarrier2(ik)-xBarrier2(ik-1);
            end
        end
       dim=length(Dx2); meanDx2=mean(Dx2(ceil(dim/6):dim));
        A2=sum(yBarrier2.*meanDx2);
        A3=(xTraf01(i-1)-xTraf02(i))*(yTraf02(i)+yTraf01(i-1)/2);
        Afe(i)=Afe(i-1)+A2+A3-A1;    % Area mezza barriera di flux
    end
    
end

    M_Fe = 2*Afe*geo.l * 1e-9 * 7800 ;                        % massa ferro appeso ai ponticelli
    
    rG=0.5*([xr,B1k]+[B2k,0]);
    rG=rG(1:nlay);
    
    F_centrifuga = M_Fe .* rG/1000 *  (nmax * pi/30)^2;
    sigma_max = 180;                                     % N/mm2 - snervamento lamierino
    
    pont = F_centrifuga/(sigma_max * geo.l);               % mm
    
    for jj=1:nlay
    if pont(jj)<pont0
        pont(jj)=0;                             % NOTA BENE: Elimino i ponticelli troppo sottili
    end
    end

