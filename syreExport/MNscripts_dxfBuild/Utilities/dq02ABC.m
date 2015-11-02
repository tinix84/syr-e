function [xA,xB,xC] = dq02ABC(theta,xD,xQ,x0)

% dq02ABC.m [v1.00.00 (30-11-2012)]
% dq0 to ABC coordinates transformation
% =========================================================================
% Syntax: [xA,xB,xC] = dq02ABC(theta,xD,xQ,x0)
% Input:
%         - theta:   angular position respect to A-phase axis [rad]
%         - xD:      d-axis component 
%         - xQ:      q-axis component 
%         - x0:      zero sequence component  
% Output:
%         - xA:      A-phase component 
%         - xB:      B-phase component  
%         - xC:      C-phase component  
% =========================================================================

nPts = length(theta);
for k = 1 : nPts
    T = [       cos(theta(k))         -sin(theta(k))  1;
        cos(theta(k)-2/3*pi)  -sin(theta(k)-2/3*pi)  1;
        cos(theta(k)-4/3*pi)  -sin(theta(k)-4/3*pi)  1];
    
    xABC = T*[xD(k);xQ(k);x0(k)];
    xA(k) = xABC(1);
    xB(k) = xABC(2);
    xC(k) = xABC(3);
end
