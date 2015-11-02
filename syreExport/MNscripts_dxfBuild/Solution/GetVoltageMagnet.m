function V = GetVoltageMagnet(magData,coilName,pID,t)

% GetVoltageMagnet.m [v1.00.00 (04-05-2013)]
% Gets instantaneous value of voltage on a coil
% =========================================================================
% Syntax: V = GetVoltageMagnet(magData,coilName,pID,t)
% Input:
%         - magData:  Magnet's data structure
%         - coilName: name of the coil
%         - pID:      Magnet's problem ID
%         - t:        time instant (-1 for magnetostatic)
%
% Output:
%         - V:        istantaneous voltage value
% =========================================================================


dh = magData.documentHandler;
mh = magData.magnetHandler;

if t<0
    % Magnetostatic
    V = invoke(dh,'getVoltageAcrossCoil',pID,coilName);
else
    % Transient -> solID = [pID t]
    invoke(mh,'processCommand','ReDim solID(1)');
    invoke(mh,'processCommand',['solID(0)=',num2str(pID)]);
    invoke(mh,'processCommand',['solID(1)=',num2str(t)]);
    
    invoke(mh,'processCommand',['CALL getDocument().getSolution().getVoltageAcrossCoil(solID,"',coilName,'", magnitude)']);
    invoke(mh,'processCommand','CALL setVariant(0,magnitude)');
    V = invoke(mh, 'getVariant', 0);
end

