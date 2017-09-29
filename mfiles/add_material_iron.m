function add_material_iron(MatName)

prompt={'Yield strength [MPa]',...
        'Density [kg*m^-3]',...
        'alpha (iron loss coeff)',...
        'beta (iron loss coeff)',...
        'kh (iron loss coeff)',...
        'ke (iron loss coeff)'};
name='New iron material';
numlines=1;
defaultanswer={'200','7800','0','0','0','0'};

answer=inputdlg(prompt,name,numlines,defaultanswer);

[filename,pathname,~]=uigetfile([cd '\.m'],'Load BH curve');
run([pathname filename]);
figure
plot(BH(:,2),BH(:,1));
if (BH(1,1)~=0 || BH(1,2)~=0)
    error('First point of BH curve must be (0,0)')
end

mat.MatName   = MatName;
mat.sigma_max = eval(answer{1});
mat.kgm3      = eval(answer{2});
mat.alpha     = eval(answer{3});
mat.beta      = eval(answer{4});
mat.kh        = eval(answer{5});
mat.ke        = eval(answer{6});
mat.BH        = BH;

load('materialLibrary\iron_material')
ind=length(MatList);
ind=ind+1;
MatList{ind}=MatName;
MatLib{ind}=mat;

button = questdlg('Save?','SELECT','Yes','No','Yes');
if isequal(button,'Yes')
    save('materialLibrary\iron_material.mat','MatList','MatLib');
    disp('material saved')
else
    disp('material not saved')
end