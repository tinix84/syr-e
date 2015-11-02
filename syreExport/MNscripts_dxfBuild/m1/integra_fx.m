% data la funzione f(x) sul dominio x, restituisce la funzione y =
% integrale(f)
% il passo su x è variabile
function y = integra_fx(x,f)

[nr numpts]=size(x);
y=0;
for kk=1:numpts
    y(kk+1)=y(kk)+f(kk);
end
y = y(2:end);

