function out=spettro(ValoreNelTempo,n,fig)

dime = max(size(ValoreNelTempo));
% non si possono richiedere più armoniche del numero di campioni-1
if (n>dime)
    n=dime-1;
end
ValoreNelTempo=reshape(ValoreNelTempo,dime,1);

a=fft(ValoreNelTempo);
Continua=a(1)/dime;
Armoniche=2*abs(a(2:dime))/dime;  % la procedura seguente serve per eliminare il valore medio plottando solo le armoniche
                                  % /dme restituisce i valori pari al
                                  % valore medio Ma sono in percentuale
                                  % rispetto alla coppia media?????????????
%Armoniche=abs(a(1:dime))/dime;
%keyboard
numarm=1:n;
% keyboard
if (fig ~= 0)
    figure(fig)
    %bar([0 numarm],[Continua;Armoniche(numarm)],0.5,'r')
    bar(numarm,Armoniche(numarm),0.5,'r')
    xlabel('0rdine Armonico'),grid,
    set(gca,'xTickLabel',numarm),
    set(gca,'xTick',1:1:length(numarm))
    
%     figure(fig+1)
%     plot(ValoreNelTempo,'b-');
%     grid on
    
end

% out=[Continua;Armoniche(numarm)];
out=Armoniche(numarm);
