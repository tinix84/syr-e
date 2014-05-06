clear res
clc
n=matlabpool('size');
if n==0
    matlabpool local 4
end
x=[6.0000    0.3333    0.3333    0.3333    0.2000    0.2000    0.2000   20.0000]';
eval_type='MO_OA';
FitnessFunction = @(x)FEMMfitnessX(x,eval_type);
n=matlabpool('size')
parfor i=1:n
    %openfemm
    res(i,:)=FitnessFunction(x);
    %closefemm
end
res
