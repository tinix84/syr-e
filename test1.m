for i=1:16
    names_o{i}=['run' num2str(i)];
end
pathname='.';
% parfor i=1:2
%     openfemm
%     tmp_fem=[pathname,'\',names_o{i} '.fem'];
%     copyfile([pathname,'\mot0.fem'],tmp_fem);
%     %h_temp=actxserver('femm.ActiveFEMM');
%     callfemm([ 'setcurrentdirectory(' , quote(pathname) , ')']);
%     opendocument(tmp_fem);
%     %FitnessFunction(x)
%     closefemm
% end
x=[6.0000    0.3333    0.3333    0.3333    0.2000    0.2000    0.2000   20.0000]';
eval_type='MO_OA';
FitnessFunction = @(x)FEMMfitness(x,eval_type);
parfor i=1:2
    %openfemm
    res(:,i)=FitnessFunction(x);
    %closefemm
end
res
