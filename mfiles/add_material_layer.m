function add_material_layer(MatName)

prompt={'Relative permeability',...
        'Density [kg*m^-3]'...
        'Remanence [T]',...
        'Coercive field [A/m]',...
        'Conductivity [S/m]'};
name='New layer material';
numlines=1;
defaultanswer={'1','0','0','0','0'};

answer=inputdlg(prompt,name,numlines,defaultanswer);

mat.MatName   = MatName;
mat.sigma     = eval(answer{1});
mat.kgm3      = eval(answer{2});

load('materialLibrary\layer_material')
ind=length(MatList);
ind=ind+1;
MatList{ind}=MatName;
MatLib{ind}=mat;

button = questdlg('Save?','SELECT','Yes','No','Yes');
if isequal(button,'Yes')
    save('materialLibrary\layer_material.mat','MatList','MatLib');
    disp('Material added to Material Library')
else
    disp('material not saved')
end