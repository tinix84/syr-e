% given the vector f(x) associated to the x-data vector x, outputs the vector y = intergal(f(x))
% x can have variable step

function y = integral_fx(x,f)

[nr numpts]=size(x);
y=0;
for kk=1:numpts
    y(kk+1)=y(kk)+f(kk);
end
y = y(2:end);

