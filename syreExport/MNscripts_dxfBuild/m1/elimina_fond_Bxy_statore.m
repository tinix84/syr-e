function Bout = elimina_fond_Bxy_statore(Bin)

[n,m] = size(Bin);
Bout = zeros(size(Bin));
x = 0:n-1;
for jj = 1:m
    temp_in = Bin(:,jj);
    temp_out = temp_in - (temp_in(1) + (temp_in(end)-temp_in(1))/n * x');
    Bout(:,jj) = temp_out;
end