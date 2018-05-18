function [FemmProblem, seginds] = addsegments_mfemm(FemmProblem, n0, n1, varargin)
% adds a segment to an mfemm FemmProblem structure.
%
% Syntax
%
% [FemmProblem, seginds] = addsegments_mfemm(FemmProblem, n0, n1)
% [FemmProblem, seginds] = addsegments_mfemm(FemmProblem, n0, n1, 'Parameter', 'Value')
% [FemmProblem, seginds] = addsegments_mfemm(FemmProblem, n0, n1, options)
%
% Description
%
% addsegments_mfemm(FemmProblem, n0, n1) adds a segment to the FemmProblem
% structure between nodes with the pairs of ids in n0 and n1.  n0 and n1
% can be matrices of node ids of the same size, with segments being added
% between each corresponding n0 and n1 id.
%
% addsegments_mfemm(FemmProblem, n0, n1, 'Parameter', 'Value') performs the
% same operation, but applies further optional parameters to the segments,
% as specified in one or more paramter-value pairs. Possible p-v pairs are:
%
%   'MaxSideLength'  - scalar value of the maximum length of triangle sides
%                      along the segment. Defaults to -1, which means no
%                      maximum length is set.
%
%   'Hidden'         - Scalar value determining the visibility in the
%                      femm post-processor. If this evaluates to true the
%                      segment will be hidden. Defaults to false.
%
%   'InGroup'        - Scalar value containing the group number of segment.
%                      Defaults to zero.
%
%   'BoundaryMarker' - String containing the name of a boundary assigned to
%                      segment. Defaults to an empty string, i.e. no
%                      boundary property.
%
% The values are applied to all segments created.
%
% addsegments_mfemm(FemmProblem, n0, n1, options) performs the
% same operation, but instead of supplying the optional arguments as p-v
% pairs they arer supplied as fields in an options structure. The fields
% must have the same names as specified in the description of the p-v pairs
% above.
%
%
% See also: newsegments_mfemm.m
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

%don't add if line is degenerate
if (n0==n1)
    return
end
if (numel(n0)>1)
    seginds=[];
    for i=1:numel(n0)
        [FemmProblem, si]=addsegments_mfemm(FemmProblem,n0(i),n1(i),varargin);
        seginds=[seginds; si];
    end
    return
end

% if (numel(varargin) == 1) && isstruct(varargin{1}) && (numel(varargin{1}) == numel(n0))
%     v=varargin{1}(1);
%     if (~isempty(varargin{1})
%    % v=varargin;
% else
%     v={};
% end

Segment.MaxSideLength = -1;
Segment.Hidden = zeros(size(n0));
Segment.InGroup = zeros(size(n0));
Segment.BoundaryMarker = '';
Segment = parseoptions(Segment, varargin);

if ~isfield(FemmProblem, 'Segments') || isempty(FemmProblem.Segments)
    seginds(1) = 1;
    
    FemmProblem.Segments = newsegment_mfemm(n0(1), n1(1),'Hidden',Segment.Hidden,...
        'InGroup',Segment.InGroup,...
        'BoundaryMarker',Segment.BoundaryMarker,...
        'MaxSideLength',Segment.MaxSideLength);
    
    linelist=FemmProblem.Segments;
    nodelist=FemmProblem.Nodes;
    dmin=norm(nodelist(n1(1)+1).Coords-nodelist(n0(1)+1).Coords)*1.e-05;
    k=numel(linelist);
    for j=0:numel(nodelist)-1
        if (j~=n0(1))&&(j~=n1(1))
            d=shortestDistance(FemmProblem,...
                nodelist(j+1).Coords(1),nodelist(j+1).Coords(2),k);
            if (norm(nodelist(j+1).Coords-nodelist(n0+1).Coords)<dmin),
                d=2.*dmin;
            end
            if (norm(nodelist(j+1).Coords-nodelist(n1+1).Coords)<dmin),
                d=2.*dmin;
            end
            if (d<dmin)
                FemmProblem.Segments(k)=[];
                FemmProblem=addsegments_mfemm(FemmProblem,n0,j,...
                    'InGroup',Segment.InGroup,'Hidden',Segment.Hidden,...
                    'MaxSideLength',-1,...
                    'BoundaryMarker',Segment.BoundaryMarker);
                FemmProblem=addsegments_mfemm(FemmProblem,j,n1,...
                    'InGroup',Segment.InGroup,'Hidden',Segment.Hidden,...
                    'MaxSideLength',-1,...
                    'BoundaryMarker',Segment.BoundaryMarker);
                j = numel(nodelist)-1;
            end
        end
    end
else
    %don't add if the line is already in the list
    for i=1:numel(FemmProblem.Segments)
        if ((FemmProblem.Segments(i).n0==n0) && (FemmProblem.Segments(i).n1==n1)), return,
        end
        if ((FemmProblem.Segments(i).n0==n1) && (FemmProblem.Segments(i).n1==n0)), return,
        end
    end
    linelist = FemmProblem.Segments;
    nodelist = FemmProblem.Nodes;
    seginds = repmat(-1, 1, numel(n0));
    % set up default segment properties
    tol=Segment.MaxSideLength;
    for i=1:numel(n0)
        %check to see if there are intersections with segments
        %If affermative, add new nodes.
        for j=1:numel(linelist)
            [xi,yi]=getIntersection(FemmProblem,n0(i),n1(i),j);
            if ~isempty(xi)
                FemmProblem=addnodes_mfemm(FemmProblem,xi,yi,'InGroup',Segment.InGroup);
            end
        end
        % Add proposed line segment
        if i>1
            seginds(i) = seginds(i-1) + 1;
        else
            seginds(i) = numel(FemmProblem.Segments) + 1;
        end
        FemmProblem.Segments(seginds(i)) = newsegment_mfemm(n0(i), n1(i),'Hidden',Segment.Hidden,...
            'InGroup',Segment.InGroup,...
            'BoundaryMarker',Segment.BoundaryMarker,...
            'MaxSideLength',Segment.MaxSideLength);
        
        % check to see if proposed line passes through other points;
        % if so, delete line and create lines that link intermediate points;
        % does this by recursive use of AddSegment;
        %if(tol==0)
        dmin=norm(nodelist(n1(i)+1).Coords-nodelist(n0(i)+1).Coords)*1.e-06;
        %else
        %    dmin=tol;
        %end
        linelist=FemmProblem.Segments;
        nodelist=FemmProblem.Nodes;
        k=numel(linelist);
        for j=0:numel(nodelist)-1
            if (j~=n0(i))&&(j~=n1(i))
                d=shortestDistance(FemmProblem,...
                    nodelist(j+1).Coords(1),nodelist(j+1).Coords(2),k);
                if (norm(nodelist(j+1).Coords-nodelist(n0+1).Coords)<dmin),
                    d=2.*dmin;
                end
                if (norm(nodelist(j+1).Coords-nodelist(n1+1).Coords)<dmin),
                    d=2.*dmin;
                end
                if (d<dmin)
                    FemmProblem.Segments(k)=[];
                    FemmProblem=addsegments_mfemm(FemmProblem,n0,j,...
                        'InGroup',Segment.InGroup,'Hidden',Segment.Hidden,...
                        'MaxSideLength',-1,...
                        'BoundaryMarker',Segment.BoundaryMarker);
                    FemmProblem=addsegments_mfemm(FemmProblem,j,n1,...
                        'InGroup',Segment.InGroup,'Hidden',Segment.Hidden,...
                        'MaxSideLength',-1,...
                        'BoundaryMarker',Segment.BoundaryMarker);
                    j = numel(nodelist)-1;
                end
            end
        end
    end
end

end


