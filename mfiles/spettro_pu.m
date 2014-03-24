function out=spettro_pu(ValoreNelTempo,n,fig)

dime = max(size(ValoreNelTempo));
% non si possono richiedere più armoniche del numero di campioni-1
if (n>dime)
    n=dime-1;
end
ValoreNelTempo=reshape(ValoreNelTempo,dime,1);

a=fft(ValoreNelTempo);
Continua=a(1)/dime;
Armoniche=2*abs(a(1:dime))/dime;
numarm=1:n;

if (fig ~= 0)
    figure(fig)
    bar(numarm,Armoniche(numarm)/Continua,2,'r')
    % bar(numarm,Armoniche(numarm),0.5,'r')
    ylabel('p.u. della continua'), xlabel('0rdine Armonico'),grid
    
    
figure(fig+1)
plot(ValoreNelTempo,'b-');
grid

end

out=[Armoniche(numarm)/Continua];
