% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function FemmProblem=draw_lines_archesX(FemmProblem,Mat,gruppo,res)
% Importa_Mat(Mat,gruppo)
[nrig,ncol] = size(Mat);
%stop=82;
for ii=1:nrig
     %if (nrig<160)&&(ii==stop)
     %    disp('X')
%         writefemmfile(['testZeroX' num2str(ii) '.fem'], FemmProblem)
%         system(['move testZeroX' num2str(ii) '.fem' ' ..'])
%         cd('..')
%         error('OK X')
   %  end
    if Mat(ii,ncol)==0
        
        outernodes = [ Mat(ii,3), Mat(ii,4); Mat(ii,1), Mat(ii,2) ];
        if(norm(outernodes(1,:)-outernodes(2,:))>eps)
            [FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2),'InGroup',gruppo);
            [FemmProblem] = addsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),...
                'MaxSideLength',res,'Hidden',0,'InGroup',gruppo,'BoundaryMarker','None');
        end

    else
        [outernodes,maxsegdeg,~,~,ang]=Disegna_ArcoX(Mat(ii,:),res);
        [FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, ...
            outernodes(:,1), outernodes(:,2),'InGroup',gruppo);
        [FemmProblem, rcsegind] = addarcsegments_mfemm(FemmProblem, nodeids(1), ...
            nodeids(2), ang, 'MaxSegDegrees', maxsegdeg, ...
            'InGroup', gruppo, 'BoundaryMarker','None');
        
    end
end
% writefemmfile('testZeroX.fem', FemmProblem)
% system(['move testZeroX.fem ..'])
% cd('..')
% error('OK X')





