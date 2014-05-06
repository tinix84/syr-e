function [centre, r] = circcentre(A, B, angle)
% calculates the centre and radius of a circle given two points and an arc
% angle between them. The position of the circle is determined by the
% order of the supplied points.
%
% Syntax
%
% [centre, r] = circcentre(A, B, angle)
%
%
    
    % get vector pointing from A to B
    AB = B - A;
    
    % find perpendicular vector to AB
    V = [ -AB(:,2), AB(:,1) ];
    
    % find mid point of AB
    M = A + AB .* 0.5;
    
    % find length of AB and divide by two to get triangle base
    b = 0.5 * magn(AB);
    
    % find triangle height
    h = b ./ tan(angle ./ 2);
    
    % find circle centre
    centre = M + h * unit(V);
    
    % find radius
    r = sqrt(h.^2 + b.^2);
    
end