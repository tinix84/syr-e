fdfq_idiq [ris1a] = interp_polin_2_(x, y, mask)
%
if ~isequal(size(x),size(y))
    error('X and Y vectors must be the same size.')
end
%
x = x(:)';
y = y(:)';
%
np = length(x);
grado = length(mask) -1;
prod = ones(grado + grado,np);
for index = 2:1:grado+grado+1;
   prod(index,:) = prod(index-1,:) .* x;
end
%
sumprod = sum(prod,2);
sumprod = flipud(sumprod);
%
m1 = prod(1:grado+1,:) .* repmat(y, grado+1,1);
sumprody = sum(m1,2);
sumprody = flipud(sumprody);
%
m2 = zeros(grado+1, grado);
index1 = 1;
for index = 1:1:grado+1;
   if (mask(1,index) ~= 0)
      m2(:,index1) = sumprod(index :(index + grado));
      index1 = index1 + 1;  
   end
end
%
index1 = 1;
for index = 1:1:grado+1;
   if (mask(1,index) ~= 0)
      m2a(index1,:) = m2(index,:);
      sumprodya(index1,:) = sumprody(index);   
      index1 = index1 +1;
   end
end
%------------------------------------------------------------------
ris1 = m2a \ sumprodya;
%
%[Q3,R3] = qr(m2a,0);
%ris1 = R3\(Q3'*sumprodya);
%
ris1a = zeros(1, grado+1);
index1 =1;
for index = 1:1:grado+1;
   if (mask(1,index) ~= 0)
      ris1a(1, index) = ris1(index1);
      index1 = index1 +1;
   end
end

