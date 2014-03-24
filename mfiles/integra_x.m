function y = integra_x(x,f)

[nr numpts]=size(x);
passo = (x(end) - x(1))/numpts;
y=0;
for kk=1:numpts
    y=y+f(kk)*passo;
end

