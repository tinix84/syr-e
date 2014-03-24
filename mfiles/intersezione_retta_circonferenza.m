%% Matteo 31/10/2011
%Intersezione retta circonferenza:
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

