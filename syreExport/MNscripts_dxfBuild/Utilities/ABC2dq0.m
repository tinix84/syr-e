function [xD,xQ,x0] = ABC2dq0(theta,xA,xB,xC)

% ABC2dq0.m [v1.00.00 (30-11-2012)]
% ABC to dq0 coordinates transformation
% =========================================================================
% Syntax: [xD,xQ,x0] = ABC2dq0(theta,xA,xB,xC)
% Input:
%         - theta:   angular position respect to A-phase axis [rad]
%         - xA:      A-phase component 
%         - xB:      B-phase component  
%         - xC:      C-phase component  
% Output:
%         - xD:      d-axis component 
%         - xQ:      q-axis component 
%         - x0:      zero sequence component  
% =========================================================================

nPts = length(theta);

for k = 1 : nPts
    T = 2/3*[ cos(theta(k))  cos(theta(k)-2/3*pi)  cos(theta(k)+2/3*pi);
             -sin(theta(k)) -sin(theta(k)-2/3*pi) -sin(theta(k)+2/3*pi);
                        1/2                   1/2                  1/2];
    
    xDQ = T*[xA(k);xB(k);xC(k)];
    xD(k) = xDQ(1);
    xQ(k) = xDQ(2);
    x0(k) = xDQ(3);
end