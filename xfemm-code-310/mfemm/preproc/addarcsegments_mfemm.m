function [FemmProblem, arcseginds] = addarcsegments_mfemm(FemmProblem, n0, n1, arcLength, varargin)
% adds arc segments to an mfemm FemmProblem structure
%
% Syntax
%
% [FemmProblem, arcseginds] = addarcsegments_mfemm(FemmProblem, n0, n1, arcLength)
% [FemmProblem, arcseginds] = addarcsegments_mfemm(FemmProblem, n0, n1, arcLength, 'Parameter', 'Value')
% [FemmProblem, arcseginds] = addarcsegments_mfemm(FemmProblem, n0, n1, arcLength, options)
%
% Description
%
% addarcsegments_mfemm(FemmProblem, n0, n1, arcLength) adds an arc segment to
% the FemmProblem structure between nodes with the pairs of ids in n0 and
% n1, sweeping out angle 'arcLength' in degrees. n0 and n1 can be matrices of
% node ids of the same size, with arc segments being added between each
% corresponding n0 and n1 id. arcLength must be a matrix or vector of the same
% size and n1 and n0.
%
% addarcsegments_mfemm(FemmProblem, n0, n1, arcLength, 'Parameter', 'Value')
% performs the same operation, but applies further optional parameters to
% the arc segments, as specified in one or more paramter-value pairs.
% Possible p-v pairs are:
%
%   'MaxSegDegrees'  - scalar value of the maximum degrees in the arc
%                      segment spanned by a single sub-segment in the
%                      discretization of the arc along the segment.
%                      Defaults to 1 degree.
%
%   'Hidden'         - Scalar value determining the visibility in the
%                      femm post-processor. If this evaluates to true the
%                      segment will be hidden. Defaults to false.
%
%   'InGroup'        - Scalar value containing the group number of segment.
%                      Defaults to zero.
%
%   'BoundaryMarker' - Either a string containing the name of a boundary
%                      assigned to segment, or an integer boundary number.
%                      If an integer, this must be the (1-base) index of an
%                      existing boundary condition in the FemmProblem
%                      Structure to be applied. If zero, no boundary
%                      property is applied. Defaults to an empty string,
%                      i.e. no boundary property.
%
% The values are applied to all segments created.
%
% addarcsegments_mfemm(FemmProblem, n0, n1, arcLength, options) performs the
% same operation, but instead of supplying the optional arguments as p-v
% pairs they arer supplied as fields in an options structure. The fields
% must have the same names as specified in the description of the p-v pairs
% above.
%
%
% See also: newarcsegment_mfemm.m
%


% Copyright 2012 Richard Crozier
% Copyright 2014 Ernesto Mininno (e.mininno@gmail.com)
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

% don't add if line is degenerate
if (n0==n1)
    return
end
% if length(n0)==1
%     nodeDistance=norm(FemmProblem.Nodes(n0+1).Coords-FemmProblem.Nodes(n1+1).Coords);
%     if (nodeDistance < 1e-6)
%         return
%     end
% end
if (numel(n0)>1)
    arcseginds=[];
    for i=1:numel(n0)
        [FemmProblem, si]=addarcsegments_mfemm(FemmProblem,n0(i),n1(i),arcLength(i),varargin);
        arcseginds=[arcseginds; si];
    end
    return
end

% set up default segment properties
options.MaxSegDegrees = ones(size(n0));
options.Hidden = zeros(size(n0));
options.InGroup = zeros(size(n0));
options.BoundaryMarker = '';
options = parseoptions(options, varargin);

if ischar(options.BoundaryMarker)
    options.BoundaryMarker = {options.BoundaryMarker};
elseif isscalar(options.BoundaryMarker)
    options.BoundaryMarker = {FemmProblem.BoundaryProps(options.BoundaryMarker).Name};
end

arcseginds = repmat(-1, 1, numel(n0));
if ~isfield(FemmProblem, 'ArcSegments') || isempty(FemmProblem.ArcSegments)
    arcseginds(1) = 1;
    FemmProblem.ArcSegments = newarcsegment_mfemm(n0, n1, arcLength, ...
        'MaxSegDegrees', options.MaxSegDegrees(1), ...
        'Hidden', options.Hidden(1), ...
        'InGroup', options.InGroup(1), ...
        'BoundaryMarker', options.BoundaryMarker{1});
else
    
    %don't add if the arc is already in the list
    for i=1:numel(FemmProblem.ArcSegments)
        if(FemmProblem.ArcSegments(i).n0==n0) && ...
                (FemmProblem.ArcSegments(i).n1==n1) && ...
                (abs(FemmProblem.ArcSegments(i).ArcLength-arcLength(1))<1.e-02),return
        end
        % arcs are ``the same'' if start and end points are the same, and if
        % the arc lengths are relatively close (but a lot farther than
        % machine precision
    end
    arcseginds(1) = numel(FemmProblem.ArcSegments) + 1;
    FemmProblem.ArcSegments(arcseginds(1)) = newarcsegment_mfemm(n0, n1, arcLength, ...
        'MaxSegDegrees', options.MaxSegDegrees(1), ...
        'Hidden', options.Hidden(1), ...
        'InGroup', options.InGroup(1), ...
        'BoundaryMarker', options.BoundaryMarker{1});
end
% add proposed arc to the arclist

arclist = FemmProblem.ArcSegments;
arc=arclist(arcseginds(1));
if isfield(FemmProblem, 'Segments')
    if ~isempty(FemmProblem.Segments)
        linelist = FemmProblem.Segments;
        for i=1:numel(linelist)
            [p,j]=getLineArcIntersection(FemmProblem, arc, linelist(i));
            if(j>0)
                for k=1:j
                    FemmProblem=addnodes_mfemm(FemmProblem,...
                        real(p(k)),imag(p(k)),'InGroup',options.InGroup(1));
                end
            end
        end
    end
end
%check if the current arc is still in the list
go=0;
for j=1:numel(arclist)
    if (FemmProblem.ArcSegments(j).n0==arc.n0) && (FemmProblem.ArcSegments(j).n1==arc.n1)
        go=1;
        break
    end
end
if go
    for i=1:numel(arclist)
        
        [p,j]=getArcArcIntersection(FemmProblem,arc,arclist(i));
        if (j>0)
            for k=1:j
                FemmProblem=addnodes_mfemm(FemmProblem,...
                    real(p(k)),imag(p(k)),'InGroup',options.InGroup(1));
            end
        end
    end
    
    % check to see if proposed arc passes through other points;
    % if so, delete arc and create arcs that link intermediate points;
    % does this by recursive use of AddArcSegment;
    %for j=1:numel(FemmProblem.ArcSegments)
    %arc=FemmProblem.ArcSegments(j);
    nodelist = FemmProblem.Nodes;
    [c,R]=getCircle(FemmProblem,arc);
    dmin=abs(R*pi*arc.ArcLength/180.)*1.e-08;
    arclist = FemmProblem.ArcSegments;
    k=numel(arclist);
    for i=1:numel(nodelist)
        if(i-1~=arc.n0)&&(i-1~=arc.n1)
            p=nodelist(i).Coords(1)+1i*nodelist(i).Coords(2);
            d=shortestDistanceFromArc(FemmProblem,arc,p);
            if d<dmin
                a0=nodelist(arc.n0+1).Coords(1)+...
                    1i*nodelist(arc.n0+1).Coords(2);
                a1=nodelist(arc.n1+1).Coords(1)+...
                    1i*nodelist(arc.n1+1).Coords(2);
                a2=nodelist(i).Coords(1)+1i*nodelist(i).Coords(2);
                FemmProblem.ArcSegments(k)=[];
                newarc=arc;
                newarc.n1=i-1;
                newarc.ArcLength=angle((a2-c)/(a0-c))*180./pi;
                FemmProblem=addarcsegments_mfemm(FemmProblem,...
                    newarc.n0,newarc.n1,newarc.ArcLength,'MaxSegDegrees', dmin, ...
                    'Hidden', options.Hidden(1), ...
                    'InGroup', options.InGroup(1), ...
                    'BoundaryMarker', options.BoundaryMarker{1});
                newarc=arc;
                newarc.n0=i-1;
                newarc.ArcLength=angle((a1-c)/(a2-c))*180./pi;
                FemmProblem=addarcsegments_mfemm(FemmProblem,...
                    newarc.n0,newarc.n1,newarc.ArcLength,'MaxSegDegrees', dmin, ...
                    'Hidden', options.Hidden(1), ...
                    'InGroup', options.InGroup(1), ...
                    'BoundaryMarker', options.BoundaryMarker{1});
            end
        end
        %  end
    end
end


end