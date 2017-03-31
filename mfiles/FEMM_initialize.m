function FEMM_initialize(geo,mat)

%% problem definition
newdocument(0);
mi_probdef(0,'millimeters','planar',1e-8,geo.l,15);

%% add iron used
mi_addmaterial(mat.Stator.MatName);

for ii=1:length(mat.Stator.BH(:,1))
    mi_addbhpoint(mat.Stator.MatName,mat.Stator.BH(ii,1),mat.Stator.BH(ii,2));
end

if ~strcmp(mat.Stator.MatName,mat.Rotor.MatName)
    mi_addmaterial(mat.Rotor.MatName);
    for ii=1:length(mat.Rotor.BH(:,1))
        mi_addbhpoint(mat.Rotor.MatName,mat.Rotor.BH(ii,1),mat.Rotor.BH(ii,2));
    end
end

if ~strcmp(mat.Stator.MatName,mat.Shaft.MatName)
    mi_addmaterial(mat.Shaft.MatName);
    for ii=1:length(mat.Shaft.BH(:,1))
        mi_addbhpoint(mat.Shaft.MatName,mat.Shaft.BH(ii,1),mat.Shaft.BH(ii,2));
    end
end

%% add conductor used
mi_addmaterial(mat.SlotCond.MatName);
mi_modifymaterial(mat.SlotCond.MatName,5,mat.SlotCond.sigma/1e6);

%% add barrier used
if isfield(mat.LayerMag,'BH')
    mi_addmaterial(mat.LayerMag.MatName);

    for ii=1:length(mat.LayerMag.BH(:,1))
        mi_addbhpoint(mat.LayerMag.MatName,mat.LayerMag.BH(ii,1),mat.LayerMag.BH(ii,2));
    end
    
else
    mi_addmaterial(mat.LayerMag.MatName,mat.LayerMag.mu,mat.LayerMag.mu,mat.LayerMag.Hc(1));
end

%% add air
mi_addmaterial('Air',1,1,0,0);