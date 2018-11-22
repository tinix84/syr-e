function remove_material(MatName,type)
% 
% remove_material(MatName,type)
% 

if strcmp(type,'Iron')
    filename='materialLibrary\iron_material.mat';
elseif strcmp(type,'Conductor')
    filename='materialLibrary\conductor_material.mat';
elseif strcmp(type,'Layer')
    filename='materialLibrary\layer_material.mat';
else
    disp('Wrong material type!')
    disp('Allowed material type:')
    disp('Iron')
    disp('Conductor')
    disp('Layer')
    error('Insert a correct material type!')
end

load(filename)

index=0;
for ii=1:length(MatList)
    if strcmp(MatName,MatList{ii})
        index=ii;
    end
end

if index==0
    disp([MatName ' not found in ' type ' material library'])
    disp('Materials present:')
    for ii=1:length(MatList)
        disp(['- ' MatList{ii}]);
    end
    error('Insert a correct material name')
end

MatListOld=MatList;
MatLibOld=MatLib;
clear MatList MatLib

ind=1;

for ii=1:length(MatListOld)
    if ii~=index
        MatList{ind}=MatListOld{ii};
        MatLib{ind}=MatLibOld{ii};
        ind=ind+1;
    end
end

save(filename,'MatList','MatLib');

disp([MatName ' deleted from ' type ' material library'])






