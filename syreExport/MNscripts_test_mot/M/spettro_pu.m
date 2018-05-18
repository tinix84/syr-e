%% VERSIONE 20 11 2011
function out=spettro_pu(ValoreNelTempo,n,fig)

dime = max(size(ValoreNelTempo));
% non si possono richiedere più armoniche del numero di campioni-1
if (n>dime)
    n=dime-1;
end
ValoreNelTempo=reshape(ValoreNelTempo,dime,1);

%figure
%plot(ValoreNelTempo,'b-');
%grid

a=fft(ValoreNelTempo);
Continua=a(1)/dime;
Armoniche=2*abs(a(2:dime))/dime;
numarm=1:n;

if (fig)
    figure(101)
    bar(numarm,Armoniche(numarm)/Continua*100,0.5,'r')
    % bar(numarm,Armoniche(numarm),0.5,'r')
    ylabel('% della continua'), xlabel('0rdine Armonico'),grid
end

out=[Armoniche(numarm)/Continua];