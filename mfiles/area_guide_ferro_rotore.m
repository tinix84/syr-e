if (ii==1)
    
%     beta = 180/pi * calc_apertura_cerchio(pi/180*alpha(ii),xr,x0);
%     r = (x0 - xr * cos(alpha(ii)*pi/180))./(cos(beta*pi/180));
%     
%     [xc,yc] = calc_intersezione_cerchi(xr-pont0-hc(ii)/2, r(ii), x0);
    
    [a1 b1 c1]=retta_tg_ad_una_circonferenza(xc(ii),yc(ii),(xc(ii)-hc(ii)/2/sqrt(2)),(yc(ii)+hc(ii)/2/sqrt(2)));
    [a2 b2 c2]=retta_per_2pti(X42(ii),Y42(ii),X32(ii+1),Y32(ii+1));
    [x5 y5]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
    [xc1,yc1]=intersezione_retta_circonferenza(xc(ii),yc(ii),xr);
    [xc2,yc2]=intersezione_retta_circonferenza(X32(ii+1),Y32(ii+1),xr);
    A5=calc_area(x5,X42(ii),(xc(ii)-hc(ii)/2/sqrt(2)),y5,Y42(ii),(yc(ii)+hc(ii)/2/sqrt(2)));% non verrà tenuto conto di questo contributo nel seguito 
    
    A4=calc_area(xc2,xc1,X32(ii+1),yc2,yc1,Y32(ii+1));
    A3=calc_area(x5,X32(ii+1),xc1,y5,Y32(ii+1),yc1);
    A2=calc_area(X42(ii),X32(ii+1),X42(ii),Y42(ii),Y32(ii+1),0);
    A1=calc_area(X32(ii+1),X32(ii+1),X42(ii),0,Y32(ii+1),0);
        
elseif(ii==nlay)
    xar=Ar; yar=Ar*tan(pi/4);
    xc2=xr*cos(pi/4);
    yc2=xr*sin(pi/4);
    [xc1,yc1]=intersezione_retta_circonferenza(X42(ii),Y42(ii),xr);
    A4=calc_area(xc2,xc1,xar,yc2,yc1,yar);
    A3=calc_area(X42(ii),xar,xc1,Y42(ii),yar,yc1);
    A2=calc_area(X42(ii),xar,X42(ii),Y42(ii),yar,0);
    A1=calc_area(xar,xar,X42(ii),0,yar,0);    

else
[xc1,yc1]=intersezione_retta_circonferenza(X42(ii),Y42(ii),xr);
[xc2,yc2]=intersezione_retta_circonferenza(X32(ii+1),Y32(ii+1),xr);
A4=calc_area(xc2,xc1,X32(ii+1),yc2,yc1,Y32(ii+1));
A3=calc_area(X42(ii),X32(ii+1),xc1,Y42(ii),Y32(ii+1),yc1);
A2=calc_area(X42(ii),X32(ii+1),X42(ii),Y42(ii),Y32(ii+1),0);
A1=calc_area(X32(ii+1),X32(ii+1),X42(ii),0,Y32(ii+1),0);    
end

Atot=A1+A2+A3+A4;