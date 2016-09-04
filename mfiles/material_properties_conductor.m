function [mat] = material_properties_conductor(MatName)


%% material list
MatList = { 'Copper'
            'Aluminium'};
mat.MatList = MatList;
mat.MatName = MatName;

if strcmp(MatName,MatList{1})
    %% Copper
    mat.sigma = 5.8e7;      % [S/m]
    mat.kgm3 = 8940;        % [kg/m3]
elseif strcmp(MatName,MatList{2})
    mat.sigma = 2.9e7;      % [S/m]
    mat.kgm3 = 2700;        % [kg/m3]
    
end
          