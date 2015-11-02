function SelectObjectFaceMagnet(magData,objectName,faceNumber,toggleFlag)

% SelectObjectFaceMagnet.m [v1.00.00 (30-11-2012)]
% Select a face of an object
% =========================================================================
% Syntax: SelectObjectFaceMagnet(magData,objectName,faceNumber,toggleFlag)
%
% Input:
%          - magData: Magnet's data structure
%          - objectName: name of the selected object
%          - faceNumber: ID number of the selected face
%          - toggleFlag: if 1 allows multiple selection (previously selected
%                        faces will not be cleared) and deselection, if 0
%                        allows only single selection and no deselection.
%
% =========================================================================


% Get View handler
vh = magData.viewHandler;

% Get Magnet's constants handler
ch = magData.magnetConstants;

% Merge object name and face number
faceName = [objectName,',Face#',num2str(faceNumber)];
% Run command
if toggleFlag==1
    invoke(vh, 'selectObject', faceName, get(ch,'infoToggleInSelection'));
else
    invoke(vh, 'selectObject', faceName, get(ch,'infoSetSelection'));
end
