function x = GetTransientSpeedAtTime(magData,ID,t,motionComponent)
% GetTransientPositionAtTime.m [v1.00.00 (02-04-2013)]
% Gets the position value for transient problems at time t.
% t could be an array to get several values.
% =========================================================================
% Syntax: x = GetTransientPositionAtTime(magData,ID,t,motionComponent)
% Input:
%          - magData: Magnet's initial data structure
%          - ID: Magnet's problem ID
%          - t: time value
%          - motionComponent: string containg the name of the motion
%                             component (e.g. "MOTION#1")
% Output:
%          - c: position value value
% =========================================================================

% Get Magnet handler 
mh = magData.magnetHandler;

% Define Val (number containing field value) and sol (array
% containing the problem ID and the time value)
invoke(mh,'processCommand','ReDim Val(0)');
invoke(mh,'processCommand','ReDim sol(1)');
invoke(mh,'processCommand',['sol(0)=',num2str(ID)]);
invoke(mh,'processCommand',['sol(1)=',num2str(t)]);

% Get the position and store it in Val
invoke(mh,'processCommand',['Val(0)=getDocument().getSolution().getSpeedOfMotionComponent(sol,',motionComponent,')']);

% Copy position value to the Variant (0,MATLAB)
invoke(mh, 'processCommand', 'call setVariant(0,Val(0),"MATLAB")');

% Get position value from the Variant (0,MATLAB)
x = invoke(mh, 'getVariant', 0, 'MATLAB');