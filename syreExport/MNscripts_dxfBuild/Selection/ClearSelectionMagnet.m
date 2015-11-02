function ClearSelectionMagnet(magData)

% ClearSelectionMagnet.m [v1.00.00 (30-11-2012)]
% Clear the current selection
% =========================================================================
% Syntax: ClearSelectionMagnet(magData)
%
% Input:
%          - magData: Magnet's data structure
%
% =========================================================================

% Get View handler
vh = magData.viewHandler;

% Run command
invoke(vh, 'unselectAll');
