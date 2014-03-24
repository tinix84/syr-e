function [front,idx] = FastParetoEstimation(x,COST) %pop,info)

pop=[x COST];
[xx,yx]=size(x);
[xC,yC]=size(COST);
info.fitSize=yC;
info.fitIdx = [yx+1:yx+info.fitSize];
[N,aux] = size(pop);
clear aux;

M = info.fitSize;
idx = ones(N,1);
fitP = zeros(1,M);
fitQ = zeros(1,M);
domFlag = 0;

for p=1:N
    if idx(p)
        fitP = pop(p,info.fitIdx);
        for q=p+1:N
            fitQ = pop(q,info.fitIdx);
            
            domFlag = Dominance(fitP,fitQ,0);
            switch domFlag
                case  1 
                    idx(q) = 0;
                case -1 
                    idx(p) = 0;
                otherwise
                        
            end
            
            if ~idx(p)
                break;
            end
        end
    end
end
                  
front = [];
for k=1:N
    if idx(k)
        front = [front;pop(k,:)];
    end
end
