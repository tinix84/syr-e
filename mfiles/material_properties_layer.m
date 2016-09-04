function [mat] = material_properties_layer(MatName)


%% material list
MatList = { 'Air'
            'Bonded-Magnet'
            'BMN-38EH'
            'BMN-42SH'
            'NdFeB 37 MGOe'};
mat.MatList = MatList;
mat.MatName = MatName;

if strcmp(MatName,MatList{1})
    %% Air
    mat.mu = 1;                 % mu_r
    mat.kgm3 = 0;              % mass density (needed for the evaluation of the radial ribs)
    mat.Br = 0;                 % Br
    mat.Hc = 0;                 % Hc
    mat.sigmaPM = 0;            % [S/m]
elseif strcmp(MatName,MatList{2})
    %% Bonded-Magnet
    mat.mu = 1;
    mat.kgm3 = 7600;
    mat.Br = 0.5;
    mat.Hc = mat.Br/(4*pi*1e-7*mat.mu);
    mat.sigmaPM = 0;
elseif strcmp(MatName,MatList{3})
    %% BMN-38EH
    mat.mu = 1.1;                       % mu_r
    mat.kgm3 = 7550;                   % [kg/m^3]
    mat.Br = 1.026;                     % [T]
    mat.Hc = mat.Br/(4*pi*1e-7*mat.mu); % [A/m]
    mat.sigmaPM = 667000;               % [S/m]
elseif strcmp(MatName,MatList{4})
    %% BMN-42SH
    mat.mu = 1.05;                      % mu_r
    mat.kgm3 = 7550;                   % [kg/m^3]
    mat.Br = 1.09;                      % [T]
    mat.Hc = mat.Br/(4*pi*1e-7*mat.mu); % [A/m]
    mat.sigmaPM = 6667000;              % [S/m]
elseif strcmp(MatName,MatList{5})
    %% NdFeB 37 MGOe
    mat.mu = 1.048;                     % mu_r
    mat.kgm3 = 7500;                   % [kg/m3]
    mat.Br = 1.251;                     % [T]
    mat.Hc = mat.Br/(4*pi*1e-7*mat.mu); % [A/m]
    mat.sigmaPM = 667000;               % [S/m]
    
    
end
          