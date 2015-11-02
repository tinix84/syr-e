F1=Flussi(:,1);F2=Flussi(:,2);F3=Flussi(:,3);

theta = (Mac.th0 * pi/180 + tempo/1000 * Cas.n * pi/30 * Mac.p);

fdq=abc2dq(F1',F2',F3',theta')
fd=fdq(1,:);fq=fdq(2,:);


thM=(theta'*180/pi-Mac.th0)/Mac.p;
% theta=(position(2:end)*Mac.p+Mac.th0)*pi/180;
figure;plot(theta,fd');

figure;plot(tempo,fq)

I=dq2abc(Id',Iq',theta')
I1=I(1,:); I2=I(2,:); I3=I(3,:);
figure;
plot(tempo',I1,'*r'); hold on;
plot(tempo',I2,'*b');
plot(tempo',I3,'*g');