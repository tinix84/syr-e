function [mat] = material_properties_layer(MatName)


%% material list
MatList = { 'Air'
            'Bonded-Magnet'
            'BMN-38EH'
            'BMN-42SH'
            'NdFeB 37 MGOe'
            'Virgin Bonded-Magnet 04T'
            'NdFeB 32 MGOe'};
mat.MatList = MatList;
mat.MatName = MatName;

if strcmp(MatName,MatList{1})
    %% Air
    mat.mu = 1;                 % mu_r
    mat.kgm3 = 0;               % mass density (needed for the evaluation of the radial ribs)
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
    mat.kgm3 = 7550;                    % [kg/m^3]
    mat.Br = 1.026;                     % [T]
    mat.Hc = mat.Br/(4*pi*1e-7*mat.mu); % [A/m]
    mat.sigmaPM = 667000;               % [S/m]
elseif strcmp(MatName,MatList{4})
    %% BMN-42SH
    mat.mu = 1.05;                      % mu_r
    mat.kgm3 = 7550;                    % [kg/m^3]
    mat.Br = 1.09;                      % [T]
    mat.Hc = mat.Br/(4*pi*1e-7*mat.mu); % [A/m]
    mat.sigmaPM = 6667000;              % [S/m]
elseif strcmp(MatName,MatList{5})
    %% NdFeB 37 MGOe
    mat.mu = 1.048;                     % mu_r
    mat.kgm3 = 7500;                    % [kg/m3]
    mat.Br = 1.251;                     % [T]
    mat.Hc = mat.Br/(4*pi*1e-7*mat.mu); % [A/m]
    mat.sigmaPM = 667000;               % [S/m]
elseif strcmp(MatName,MatList{6})
    %% Virgin Bonded-Magnet 04T
    mat.kgm3 = 7600;                    % [kg/m3]
    mat.Br = 0;
    mat.Hc = 0;
    mat.BH = [0.000000 	 0.000000
              0.013440 	 92309.298661
              0.052400 	 135183.807845
              0.103600 	 155435.990940
              0.162400 	 182370.317119
              0.230400 	 211134.351303
              0.289200 	 238917.301706
              0.328800 	 274369.660419
              0.349200 	 310564.675711
              0.362800 	 348319.181512
              0.372000 	 379002.113134
              0.380400 	 419361.656876
              0.388000 	 465809.924890
              0.393600 	 540268.090679
              0.400000 	 710000.000000
              0.400000 	 879731.909321];  %this is the characteristic J,H
    mat.BH = [mat.BH(:,1)+4*pi*1e-7*mat.BH(:,2) mat.BH(:,2)];
    mat.Bnom=0.4;  %[T]
    mat.Hci=284000; %[A/m]
elseif strcmp(MatName,MatList{7})
        %% NdFeB 32 MGOe
    mat.mu = 1.045;                     % mu_r
    mat.kgm3 = 7500;                    % [kg/m3]
    mat.Br = 1.16;                     % [T]
    mat.Hc = mat.Br/(4*pi*1e-7*mat.mu); % [A/m]
    mat.sigmaPM = 667000;               % [S/m]
    
end
          