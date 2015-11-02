function Pj = GetOhmicLossesMagnet(magData,ID,t,compName)

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
invoke(mh,'processCommand','ReDim ohmLoss(0)');

if t<0
    % Magnetostatic ID = pID
    invoke(mh,'processCommand','ReDim ID(0)');
    invoke(mh,'processCommand',['ID(0)=',num2str(ID)]);
else
    % Transient -> ID = [pID t]
    invoke(mh,'processCommand','ReDim ID(1)');
    invoke(mh,'processCommand',['ID(0)=',num2str(ID)]);
    invoke(mh,'processCommand',['ID(1)=',num2str(t)]);
end

invoke(mh,'processCommand',['ohmLoss(0)=getDocument().getSolution().getOhmicLossInConductor(ID,"',compName,'")']);
invoke(mh,'processCommand','CALL setVariant(0,ohmLoss)');
Pj = invoke(mh,'getVariant',0);
Pj = Pj{1};