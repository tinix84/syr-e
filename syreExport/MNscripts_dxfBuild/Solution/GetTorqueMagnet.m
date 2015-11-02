function T = GetTorqueMagnet(magData,center,pID,t)

% GetTorqueMagnet.m [v1.01.00 (04-04-2013)]
% Gets torque from each body
% ==========================================================================
% Syntax: T = GetTorqueMagnet(magData,center)
% Input:
%         - magData:  Magnet's data structure
%         - center:   center from which the torque is computed
%
% Output:
%         - T:        torque values vector (each components is referred to a
%                     body)
% =========================================================================
% Change Log
% 04-04-2013: added transient+motion capability

%dh = magData.documentHandler;
mh = magData.magnetHandler;
sh = magData.solutionHandler;

% Number of bodies
nBodies = invoke(sh,'getNumberOfBodies',1);
T = zeros(nBodies,1);
% Center string
centerString = ['Array(',num2str(center(1)),',',num2str(center(2)),',',num2str(center(3)),')'];

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
    invoke(mh,'processCommand',['CALL getDocument().getSolution().getTorqueOnBody(solID,',num2str(k),',',centerString,',torque_x,torque_y,torque_z)']);
    invoke(mh,'processCommand','CALL setVariant(0,torque_z)');
    T(k) = invoke(mh, 'getVariant', 0);
end

