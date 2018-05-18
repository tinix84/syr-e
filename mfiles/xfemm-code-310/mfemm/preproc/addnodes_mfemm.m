function [FemmProblem, nodeinds, nodeids] = addnodes_mfemm(FemmProblem, x, y, varargin)
% adds one or more nodes to an mfemm FemmProblem structure at the specified
% locations
%
% Syntax
%
% [FemmProblem, nodeinds, nodeids] = addnodes_mfemm(FemmProblem, x, y, varargin)
%
% Description
%
%
% Input
%
% FemmProblem - A FemmProblem structure as created by newproblem_mfemm.m
%
% x - A matrix of x (or r for axisymmetric problems)coordinates of the
%   nodes to be added to the problem structure
%
% y - A matrix of y (or z for axisymmetric problems) coordinates of the
%   nodes to be added to the problem structure. Must be the same size as x.
%
%

% Copyright 2012 Richard Crozier
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

distance=1e-8;
nodeinds = repmat(-1, 1, numel(x));
if ~isfield(FemmProblem, 'Nodes') || isempty(FemmProblem.Nodes)
    
    FemmProblem.Nodes = newnode_mfemm(x(1), y(1), varargin{:});
    nodeinds(1) = 1;
    
end
%nodeinds(1) = numel(FemmProblem.Nodes) + 1;
%FemmProblem.Nodes(nodeinds(1)) = newnode_mfemm(x(1), y(1), varargin{:});

for i=1:numel(x)
    [id, xy]=findnode_mfemm(FemmProblem,[x(i) y(i)]);
    % test to see if ``too close'' to existing node...
    if(norm(xy-[x(i) y(i)])>distance)
        nodeinds(i) = numel(FemmProblem.Nodes) + 1;
        FemmProblem.Nodes(nodeinds(i)) = newnode_mfemm(x(i), y(i), varargin{:});
    else
        nodeinds(i) = id+1;
        nodeids=nodeinds-1;
        continue
    end
    % can't put a node on top of a block label; do same sort of test.
    if isfield(FemmProblem,'BlockLabels')
        for j=1:numel(FemmProblem.BlockLabels)
            d=norm([x(i) y(i)]-FemmProblem.BlockLabels(j).Coords);
            if d<distance
                x(i)=[];
                y(i)=[];
            end
        end
    end
    if isfield(FemmProblem,'Segments')
        linelist=FemmProblem.Segments;
        nodelist=FemmProblem.Nodes;
        %test to see if node is on an existing line; if so, break into two lines;
        k=numel(linelist);
        for j=1:k
            if (abs(shortestDistance(FemmProblem,x(i),y(i),j)<distance))
                segm=linelist(j);
                FemmProblem.Segments(j).n1=numel(nodelist)-1;
                segm.n0=numel(nodelist)-1;
                FemmProblem=addsegments_mfemm(FemmProblem,segm.n0,segm.n1);
            end
        end
    end
    if isfield(FemmProblem,'ArcSegments')
        arclist=FemmProblem.ArcSegments;
        nodelist=FemmProblem.Nodes;
        k=numel(arclist);
        for j=1:k
            if(shortestDistanceFromArc(FemmProblem,arclist(j),x(i)+1i*y(i))<distance)
                xx=nodelist(arclist(j).n0+1).Coords;
                a0=xx(1)+1i*xx(2);
                xx=nodelist(arclist(j).n1+1).Coords;
                a1=xx(1)+1i*xx(2);
                a2=x(i)+1i*y(i);
                [c,~]=getCircle(FemmProblem,arclist(j));
                arc=arclist(j);
                FemmProblem.ArcSegments(j).n1=numel(nodelist)-1;
                FemmProblem.ArcSegments(j).ArcLength=angle((a2-c)/(a0-c))*180./pi;
                arc.n0=numel(nodelist)-1;
                arc.ArcLength=angle((a1-c)/(a2-c))*180./pi;
                FemmProblem=addarcsegments_mfemm(FemmProblem,...
                    arc.n0,arc.n1,arc.ArcLength);
            end
        end
    end
    nodeids = nodeinds - 1;
end