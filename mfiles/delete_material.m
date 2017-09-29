function delete_material(MatName)

load('materialLibrary\iron_material.mat')

ind=0;

for ii=1:length(MatList)
    if strcmp(MatList{ii},MatName)
        ind=ii;
    end
end
type=1;

if ind==0
    load('materialLibrary\layer_material.mat')
    for ii=1:length(MatList)
        if strcmp(MatList{ii},MatName)
            ind=ii;
        end
    end
    type=2;
end

if ind==0
    load('materialLibrary\conductor_material.mat')
    for ii=1:length(MatList)
        if strcmp(MatList{ii},MatName)
            ind=ii;
        end
    end
    type=3;
end

if ind==0
    disp('Material not saved in the library')
else
    
    indOld=0;
    
    for ii=1:length(MatList)
        if strcmp(MatList{ii},MatName)
            indOld=ii;
        end
    end
    if indOld==0
        mat.MatName=MatName;
        mat.MatList=MatList;
    else
        mat=MatLib{indOld};
    end
    
    MatListOld=MatList;
    MatLibOld=MatLib;
    clear MatList MatLib
    
    ind=1;
    
    if indOld==1
        for ii=2:length(MatListOld)
            MatList{ind}=MatListOld{ii};
            MatLib{ind}=MatLibOld{ii};
            ind=ind+1;
        end
    elseif indOld==length(MatListOld)
        for ii=1:length(MatListOld)-1
            MatList{ind}=MatListOld{ii};
            MatLib{ind}=MatLibOld{ii};
            ind=ind+1;
        end
    else
        for ii=1:indOld-1
            MatList{ind}=MatListOld{ii};
            MatLib{ind}=MatLibOld{ii};
            ind=ind+1;
        end
        for ii=indOld+1:length(MatListOld)
            MatList{ind}=MatListOld{ii};
            MatLib{ind}=MatLibOld{ii};
            ind=ind+1;
        end
    end
    
    switch type
        case 1
            save('materialLibrary\iron_material.mat')
            disp('Deleted iron material')
        case 2
            save('materialLibrary\layer_material.mat')
            disp('Deleted iron material')
        case 3
            save('materialLibrary\conductor_material.mat')
            disp('Deleted iron material')
    end
    
    
    
end

