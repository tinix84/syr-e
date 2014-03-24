% definizione della matrice dei nodi master
% statore;    % carica parametri statore
[x3,y3]=pol2cart(s.AS1*pi/180*acs,s.RSI);
[x4,y4]=pol2cart(s.AS1*pi/180,s.RSI);
[x7,y7]=pol2cart(s.AS1*pi/180,s.RS3);
[x10,y10]=pol2cart(s.AS1*pi/180,s.RS1);
[x13,y13]=pol2cart(s.AS1*pi/180,s.RS2);
[x16,y16]=pol2cart(s.AS1*pi/180,s.RSE);
[x19,y19]=pol2cart(s.AS1*pi/180,s.RS4);
[x20,y20]=pol2cart(s.AS1*pi/180,s.RS5);

mast=[];
mast=[0 0       %1
   s.RSI 0
   x3 y3
   x4 y4
   s.RSI+s.DS2 0    %5
   s.RSI+s.DS2 s.DS1
   x7 y7    %7
   s.RS1 0
   s.RS1 s.RC1
   x10 y10  %10
   s.RS2 0
   s.RS2 s.RC2
   x13 y13
   s.RS5 0
   s.RSE 0
   x16 y16
   s.RS4 0
   s.RS4 s.RC3  %18
   x19 y19
   x20 y20
   0 0
   0 0
   s.RS5 s.DS3];    %23
nummast=length(mast);





