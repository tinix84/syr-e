n=newnode_mfemm(x(1),y(1),varargin);

for i=1:numel(x)
    [id, xy]=findnode_mfemm(FemmProblem,[x(i) y(i)]);
    if(norm(xy-[x(i) y(i)])>eps)||(n.InGroup~=FemmProblem.Nodes(id+1).InGroup)
        
end