function SelectObjectMagnet(magData,objectName,toggleFlag)

% SelectObjectMagnet.m [v1.00.00 (30-11-2012)]
% Selects an object
% =========================================================================
% Syntax: SelectObjectMagnet(magData,objectName)
%
% Input:
%          - magData:    Magnet's data structure
%          - objectName: name of the object to be selected
%          - toggleFlag: if 1 allows deselection, if 0 allows selection
%                        only
% =========================================================================

% Get View handler
vh = magData.viewHandler;

% Get Magnet's constants handler
ch = magData.magnetConstants;

% Run command
if toggleFlag==1
    invoke(vh, 'selectObject', objectName, get(ch,'infoToggleInSelection'));
else
    invoke(vh, 'selectObject', objectName, get(ch,'infoSetSelection'));
end