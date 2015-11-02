function [P_hyst,P_eddy] = GetIronLossesMagnet(magData,ID,compName)

% GetOhmicLossesMagnet.m [v1.01.00 (05-04-2013)]
% Gets istantaneous value of ohmic losses on a component
% =========================================================================
% Syntax: Pj = GetOhmicLossesMagnet(magData,ID,t,compName)
% Input:
%         - magData:    Magnet's data structure
%         - pID:        Magnet's problem ID
%         - t:          time instant (-1 for magnetostatic)
%         - compName:   string containing the name of the component on which
%                       losses are computed
% Output:
%         - Pj:         istantaneous losses value [W]
% =========================================================================

mh = magData.magnetHandler;

invoke(mh,'processCommand','ReDim components(0)');
invoke(mh,'processCommand',['components(0)="',compName,'"']);

invoke(mh,'processCommand',['CALL getDocument().getSolution().getIronLossInComponent(',num2str(ID),',components,Losses)']);
invoke(mh,'processCommand','CALL setVariant(0,Losses)');
P_iron = invoke(mh,'getVariant',0);
P_hyst = P_iron{1};
P_eddy = P_iron{2};

