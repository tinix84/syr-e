% mean barrier value:
function [xBkMean,yBkMean]=valore_medio_di_barriera(xxBk,yyBk,yTraf)

xAux=zeros(1,50); yAux=zeros(1,50);
yBkMean=zeros(1,size(yyBk,1));
xBkMean=zeros(1,size(yyBk,1));
for k=1:size(yyBk,1)
%     keyboard
    pos=find(yyBk(k,:)>0 & yyBk(k,:)<=yTraf(k));
    if (length(pos)==1|| isempty(pos))
        pos=find(yyBk(k,:)>=0); 
        yAux=linspace(yyBk(k,pos(1)),yTraf(k),100);
        xAux=interp1(yyBk(k,pos),xxBk(k,pos),yAux);
        yBkMean(k)=mean(yAux);
        xBkMean(k)=interp1(yAux,xAux,mean(yAux));
    else
        yAux=linspace(yyBk(k,pos(1)),yyBk(k,pos(end)),100);
        xAux=interp1(yyBk(k,pos),xxBk(k,pos),yAux);
        yBkMean(k)=mean(yAux);
        xBkMean(k)=interp1(yAux,xAux,mean(yAux));
    end
    
end

end