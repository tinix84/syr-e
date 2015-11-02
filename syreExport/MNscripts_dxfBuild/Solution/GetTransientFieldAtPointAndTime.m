function F = GetTransientFieldAtPointAndTime(magData,ID,t,p,fieldType)
% GetTransientFieldAtPointAndTime.m [v1.00.00 (02-04-2013)]
% Gets the field value for transient problems on point p=(x,y) at time t.
% p could be a matrix [x1,y1;x2,y2;...] to sample several points.
% =========================================================================
% Syntax: F = GetTransientFieldAtPointAndTime(magData,ID,t,x,y,fieldType)
% Input:
%          - magData: Magnet's initial data structure
%          - ID: Magnet's problem ID
%          - t: time value
%          - p: point to be sampled
%          - field type: string containg the field to be captured
%                        (e.g. "|B| smoothed")
% Output:
%          - F: Field value
% =========================================================================

% Get Magnet's handler
mn6 = magData.magnetHandler;

% Define FieldVals (number containing field value) and solution (array
% containing the problem ID and the time value)
invoke(mn6, 'processCommand','ReDim FieldVals(0)');
invoke(mn6, 'processCommand','ReDim solution(1)');
invoke(mn6, 'processCommand',['solution(0)=',num2str(ID)]);
invoke(mn6, 'processCommand',['solution(1)=',num2str(t)]);

% Get the handler to the field (stored in internal variable Field)
invoke(mn6, 'processCommand',['Set Field=getDocument().getSolution().getTransientField(solution,',fieldType,')']);


[nPts,aux] = size(p);
F = zeros(nPts,1);
for k = 1 : nPts
    x = p(k,1);
    y = p(k,2);
    % Get the field at (x,y) and store it in FieldVal
    command = ['Call Field.getFieldAtPoint(',num2str(x),',',num2str(y),',0 ,FieldVals)'];
    invoke(mn6, 'processCommand', command);
    % Transfer field value to the Variant (0,MATLAB)
    invoke(mn6, 'processCommand', 'call setVariant(0,FieldVals(0),"MATLAB")');
    % Get the field value from the Variant (0,MATLAB)
    F(k) = invoke(mn6, 'getVariant', 0, 'MATLAB');
end

