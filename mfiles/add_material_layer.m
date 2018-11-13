function add_material_layer(MatName)

% prompt={'Relative permeability',...
%         'Density [kg*m^-3]'...
%         'Remanence [T]',...
%         'Coercive field [A/m]',...
%         'Conductivity [S/m]',...
%         'Reference temperature [�C]',...
%         'Br=f(temperature) [T]',...
%         'Br=f(temperature) [T]'};
prompt={'Relative permeability',...
        'Density [kg*m^-3]'...
        'Remanence [T]',...
        'Conductivity [S/m]',...
        'Reference temperature [�C]',...
        'Br=f(temperature) [T]',...
        'Bd=f(temperature) [T]'};    
name='New layer material';
numlines=1;
defaultanswer={'1','0','0','0','[20]','[0]','[0]'};
% defaultanswer={'1','0','0','0','0','[20]','[0]','[0]'};

answer=inputdlg(prompt,name,numlines,defaultanswer);

mat.MatName   = MatName;
mat.mu        = eval(answer{1});
mat.kgm3      = eval(answer{2});
mat.Br        = eval(answer{3});
%mat.Hc        = eval(answer{3});
mat.sigmaPM   = eval(answer{4});
mat.temp.temp = eval(answer{5});
mat.temp.Br   = eval(answer{6});
mat.temp.Bd   = eval(answer{7});
% mat.Hc        = eval(answer{4});
% mat.sigmaPM   = eval(answer{5});
% mat.temp.temp = eval(answer{6});
% mat.temp.Br   = eval(answer{7});
% mat.temp.Bd   = eval(answer{8});

mat.Hc = mat.Br/(mat.mu*(4e-7*pi));

load('materialLibrary\layer_material')
ind=length(MatList);
ind=ind+1;
MatList{ind}=MatName;
MatLib{ind}=mat;

figure()
figSetting()
xlabel('$\theta$ [$^\circ$C]')
ylabel('[T]')
plot(mat.temp.temp,mat.temp.Br,'-go','DisplayName','$B_r$');
plot(mat.temp.temp,mat.temp.Bd,'-ro','DisplayName','$B_d$');
legend('show')

button = questdlg('Save?','SELECT','Yes','No','Yes');
if isequal(button,'Yes')
    save('materialLibrary\layer_material.mat','MatList','MatLib');
    disp('Material added to Material Library')
else
    disp('material not saved')
end