
n=matlabpool('size')
parfor i=1:n
    currentDir{i}=pwd;
    name{i}=getComputerName;
    cList{i}=ls;
end