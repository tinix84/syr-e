function [x, y] = arcpoints(A, B, angle, maxdeg)
% get the points on an arc for plotting

    % convert the angles to radians
    angle = deg2rad(angle);
    maxdeg = deg2rad(maxdeg);

    % get centre and radius of circles
    [centre, r] = circcentre(A, B, angle);
    
    % get starting angle of arcs
    tempA = A - centre;
    
    [starttheta, rho] = cart2pol(tempA(:,1), tempA(:,2));
    
    % get arc points
    npnts = max(3, ceil(angle ./ maxdeg));
    
    pnts = linspace(starttheta, starttheta + angle, npnts);
    
    [x, y] = pol2cart(pnts, repmat(rho, size(pnts)));
    
    x = x + centre(:,1);
    y = y + centre(:,2);
    
    
end