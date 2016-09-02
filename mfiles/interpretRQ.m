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

function [geo,gamma,mat] = interpretRQ(RQ,geo,mat)

% INTERPRETRQ interprets the string of the MOOA inputs (RQ)
% and returns dalpha, hc_pu, dx, gamma
% dalpha per unit

if not(isempty(geo.RQnames))
    first_index = 1;
    if strcmp(geo.RQnames{first_index},'dalpha')
        last_index = first_index + geo.nlay-1;
        geo.dalpha_pu = RQ(first_index:last_index);
    else
        last_index = 0;
    end
    
    % the sum of pu angles is rescaled to one
    % alpha1 is not rescaled: only the angles from the second barrier onwards
    if sum(geo.dalpha_pu) > (1-0.05) % 0.05 is the angular space guaranteed for the spider
        geo.dalpha_pu(2:end) = geo.dalpha_pu(2:end)/sum(geo.dalpha_pu(2:end))*(1-0.05-geo.dalpha_pu(1));
    end
    
    first_index = last_index + 1;
    if length(geo.RQnames)>(last_index+1)
        if strcmp(geo.RQnames{first_index},'hc')
            % hc per unit
            last_index = first_index + geo.nlay - 1;
            size(RQ);
            getComputerName();
            geo.hc_pu = RQ(first_index:last_index);
            % debug
            geo.hc_pu = sort(geo.hc_pu);
        end
    end
    
    first_index = last_index + 1;
    if length(geo.RQnames)>(last_index+1)
        if strcmp(geo.RQnames{first_index},'dx')
            % dx per unit
            last_index = first_index + geo.nlay - 1;
            geo.dx = RQ(first_index:last_index);
        end
    end
    
    first_index = last_index + 1;
    if length(geo.RQnames)>(last_index+1)
        if strcmp(geo.RQnames{first_index},'Br')
            % Br [T]
            last_index = first_index + geo.nlay - 1;
            mat.LayerMag.Br = RQ(first_index:last_index);
        end
    end
    
    if length(geo.RQnames)>(last_index)
        for k = last_index+1:length(geo.RQnames)
            eval(['geo.' geo.RQnames{k} ' = ' num2str(RQ(k)) ';'])
        end
    end
    
end

geo.dalpha = geo.dalpha_pu*(90/geo.p);  % [mec degrees]

% current phase angle
if not(isempty(RQ))
    gamma = RQ(end);
else
    gamma = 45;
end

if strcmp(geo.RotType,'SPM')&&ismember('hc',geo.RQnames)
    geo.lm = geo.hc_pu;
end