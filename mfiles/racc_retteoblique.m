%Function "racc_retteoblique" Raccordo tra due rette che formano un angolo
%qualsiasi, con un raggio di circonferenza di raggio racc_pont dato
%INPUT:x1,y1,x2,y2, coordinate punti appartenenti alla retta lato
%barriera,punto P appartenente alla retta stessa ma traslata per avere
%raggio di raccordo voluto ed infine r_racc, raggio raccordo arco di circonferenza
%OUTPUT:xc,yc, coordinate centro circonferenza di raccordo e A(xA,yA),B(xB,yB) punti 
%raccordo con le rette date a partire dal centro C

function [xc,yc,rc,xA,yA,xB,yB]=racc_retteoblique(x1,y1,x2,y2,xp,yp,angle,r_racc)
rc=r_racc; %raggio del raccordo fissato
[a,b,c]=retta_per_2pti(x1,y1,x2,y2); %retta lato barriera di flusso
[mt,qt]=retta_abc2mq(a,b,c);
%Definisco due rette parallele a quelle di partenza ma traslate di rc e
%ricavo il centro della circonferenza di raccordo
xc=(mt*xp+rc)/mt; %coordinate centro arco di circonferenza di raccordo
yc=yp+rc;
%Trovo i punti di raccordo A e B a partire dal centro C, tracciando le
%perpendicolari alla rette date condotte a partire da C
xA=xp+((rc/sin(angle))*cos(angle));
yA=yp;

mt_orto=-1/mt; %coeff.angolare perpendicolare
xB=(yc-mt_orto*xc-qt)/(mt-mt_orto);
yB=mt*xB+qt;

end