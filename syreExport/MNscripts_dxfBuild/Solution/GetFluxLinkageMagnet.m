function lambda = GetFluxLinkageMagnet(magData,coilName,pID,t)

% GetFluxLinkageMagnet.m [v1.01.00 (30-11-2012)]
% Gets flux linkage through a coil
% =========================================================================
% Syntax: lambda = GetFluxLinkageMagnet(magData,coilName)
% Input:
%         - magData:  Magnet's data structure
%         - coilName: name of the coil
%         - pID:      Magnet's problem ID
%         - t:        time instant (-1 for magnetostatic)
%
% Output:
%         - lambda:   flux linkage value
% =========================================================================


dh = magData.documentHandler;
mh = magData.magnetHandler;

if t<0
    % Magnetostatic
    lambda = invoke(dh,'getFluxLinkageThroughCoil',pID,coilName);
else
    % Transient -> solID = [pID t]
    invoke(mh,'processCommand','ReDim solID(1)');
    invoke(mh,'processCommand',['solID(0)=',num2str(pID)]);
    invoke(mh,'processCommand',['solID(1)=',num2str(t)]);
    
    invoke(mh,'processCommand',['CALL getDocument().getSolution().getFluxLinkageThroughCoil(solID,"',coilName,'", magnitude)']);
    invoke(mh,'processCommand','CALL setVariant(0,magnitude)');
    lambda = invoke(mh, 'getVariant', 0);
end
