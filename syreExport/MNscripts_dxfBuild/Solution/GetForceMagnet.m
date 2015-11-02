function F = GetForceMagnet(magData,pID,t,fDir)

% GetForceMagnet.m [v1.00.00 05-04-2013)]
% Gets force respect to x, y or z direction
% ==========================================================================
% Syntax: F = GetForceMagnet(magData,pID,t,fDir)
% Input:
%         - magData:  Magnet's data structure
%         - pID:      Magnet's problem ID
%         - t:        time instant (-1 for magnetostatic)
%         - fDir:     string containing the force direction ('x','y' or 'z')
%
% Output:
%         - F:        vector of forces (each component is related to a body)
% ==========================================================================


mh = magData.magnetHandler;
sh = magData.solutionHandler;

% Number of bodies
nBodies = invoke(sh,'getNumberOfBodies',1);
F = zeros(nBodies,1);

if t<0
    % Magnetostatic -> solID = pID
    invoke(mh,'processCommand','ReDim solID(0)');
    invoke(mh,'processCommand',['solID(0)=',num2str(pID)]);
else
    % Transient -> solID = [pID t]
    invoke(mh,'processCommand','ReDim solID(1)');
    invoke(mh,'processCommand',['solID(0)=',num2str(pID)]);
    invoke(mh,'processCommand',['solID(1)=',num2str(t)]);
end

for k = 1 : nBodies
    invoke(mh,'processCommand',['CALL getDocument().getSolution().getForceOnBody(solID,',num2str(k),',force_x,force_y,force_z)']);
    
    if isequal(fDir,'x')
        invoke(mh,'processCommand','CALL setVariant(0,force_x)');
    elseif isequal(fDir,'y')
        invoke(mh,'processCommand','CALL setVariant(0,force_y)');
    elseif isequal(fDir,'z')
        invoke(mh,'processCommand','CALL setVariant(0,force_z)');
    end
    F(k) = invoke(mh, 'getVariant', 0);
end

