function y = integra(x)

[nr numpts]=size(x);
passo=2*pi/numpts;
y=0;
for kk=1:numpts
    y=y+x(kk)*passo;
end

