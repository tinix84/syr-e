function add_material_conductor(MatName)

prompt={'Conductivity [S/m]',...
        'Density [kg*m^-3]'};
name='New iron material';
numlines=1;
defaultanswer={'580000000','8900'};

answer=inputdlg(prompt,name,numlines,defaultanswer);

mat.MatName   = MatName;
mat.sigma     = eval(answer{1});
mat.kgm3      = eval(answer{2});

load('materialLibrary\conductor_material')
ind=length(MatList);
ind=ind+1;
MatList{ind}=MatName;
MatLib{ind}=mat;

button = questdlg('Save?','SELECT','Yes','No','Yes');
if isequal(button,'Yes')
    save('materialLibrary\conductor_material.mat','MatList','MatLib');
    disp('Material added to Material Library')
else
    disp('material not saved')
end