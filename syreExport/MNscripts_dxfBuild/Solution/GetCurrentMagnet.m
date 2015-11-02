function I = GetCurrentMagnet(magData,coilName,pID,t)

% GetCurrentMagnet.m [v1.00.00 (04-05-2013)]
% Gets instantaneous value of current flowing in a coil
% =========================================================================
% Syntax: I = GetCurrentMagnet(magData,coilName,pID,t)
% Input:
%         - magData:  Magnet's data structure
%         - coilName: name of the coil
%         - pID:      Magnet's problem ID
%         - t:        time instant (-1 for magnetostatic)
%
% Output:
%         - I:        istantaneous current value
% =========================================================================

dh = magData.documentHandler;
mh = magData.magnetHandler;

if t<0
    % Magnetostatic
    I = invoke(dh,'getCurrentThroughCoil',pID,coilName);
else
    % Transient -> solID = [pID t]
    invoke(mh,'processCommand','ReDim solID(1)');
    invoke(mh,'processCommand',['solID(0)=',num2str(pID)]);
    invoke(mh,'processCommand',['solID(1)=',num2str(t)]);
    
    invoke(mh,'processCommand',['CALL getDocument().getSolution().getCurrentThroughCoil(solID,"',coilName,'", magnitude)']);
    invoke(mh,'processCommand','CALL setVariant(0,magnitude)');
    I = invoke(mh, 'getVariant', 0);
end

