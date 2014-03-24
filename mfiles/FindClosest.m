function closerIdx = FindClosest(set1,set2)

[size1,tmp] = size(set1);
[size2,tmp] = size(set2);

distMatrix = zeros(size1,size2);
closerIdx = zeros(size1,1);
for i = 1 : size1
    v1 = set1(i,:);
    
    for j = 1 : size2
        
        v2 = set2(j,:);
        
        distMatrix(i,j) = euclidean(v1,v2);
        
    end
    
    [minVal,minIdx] = min(distMatrix(i,:));
    closerIdx(i) = minIdx;
    
end



function dist = euclidean(v1,v2)

d = 0;
for k = 1 : length(v1)
    d = d + (v1(k) - v2(k)).^2;
end
dist = sqrt(d);