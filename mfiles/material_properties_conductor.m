function [mat] = material_properties_conductor(MatName)

load('materialLibrary\conductor_material.mat')

ind=0;

for ii=1:length(MatList)
    if strcmp(MatList{ii},MatName)
        ind=ii;
    end
end
if ind==0
    mat.MatName=MatName;
    mat.MatList=MatList;
else
    mat=MatLib{ind};
end