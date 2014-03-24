function [x,y]=interp_flux_barrier(xxBk,yyBk,yEnd)

x=zeros(size(yyBk,1),50); y=zeros(size(yyBk,1),50);
xAux=zeros(1,50); yAux=zeros(1,50);

for k=1:size(yyBk,1)
    pos=find(yyBk(k,:)<=yEnd(k)*20);
    yAux=linspace(0,yEnd(k),50);
    xAux=interp1(yyBk(k,pos),xxBk(k,pos),yAux);
    x(k,:)=xAux;
    y(k,:)=yAux;
end
end