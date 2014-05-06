function c=bendcountour(c,angle, anglestep)
[~,col]=size(c);
contour=c(:,1)+1i*c(:,2);
if angle==0
    return
end
if anglestep==0
    anglestep=1;
end
% check to see if there are at least enough
% points to have made one line;
k= numel(contour);
if (k<2)
    return;
end
%restrict the angle of the contour to 180 degrees;
if ((angle<-180.) || (angle>180.))
    return
end
n=ceil(abs(angle/anglestep));
tta=angle*pi/180.;
dtta=tta/n;
% pop last point off of the contour;
a1=contour(k);
contour(k)=[];
a0=contour(k-1);
% compute location of arc center;
% and radius of the circle that the
% arc lives on.
d=abs(a1-a0);
R=d/(2.*sin(abs(tta/2.)));
if(tta>0)
    c=a0 + (R/d)*(a1-a0)*exp(1i*(pi-tta)/2.);
else
    c=a0+(R/d)*(a1-a0)*exp(-1i*(pi+tta)/2.);
end

% add the points on the contour
for k=1:n 
    contour=[contour; (c+(a0-c)*exp(k*1i*dtta))];
end
c=[real(contour) imag(contour)];