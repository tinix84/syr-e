%% Matteo 31/10/2011
%Intersezione retta circonferenza:
%la circonferenza in questione ha centro nell'origine e raggio di valore xr
%la retta passa per il punto x0, y0 e ha coef angolare pari a 45°.
%
function [x y]=intersezione_retta_circonferenza(x0,y0,xr)

[v]=roots([1, 1*(y0-x0), 1/2*(y0^2+x0^2-2*y0*x0-xr^2)]);

t=1;
for i=1:length(v)
    if v(i)>=0
        w(t)=v(i);
    else
        
    end
    t=t+1;
end
x=w;
y=y0+x-x0;

