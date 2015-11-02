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

function BarName = defineBlockNames(temp,geo)

% names of rotor barriers (not used in FEMM, needed for
% compatibility with other finite-element codes such as Infolytica/Magnet)

ps = geo.ps;
nlay = geo.nlay;

switch geo.RotType
    
    case 'SPM'
        
        RotName   = {'rotor'};
        ShaftName = {'shaft'};
        PMName    = {'permanent magnet'};
        kkk=1;
        BarName = [];
        for kk = 1:ps
            BarName{kkk} = {['PM_',num2str(kk)]};
            kkk=kkk+1;
        end
        
    otherwise
        
        YpontRadSx = temp.YpontRadSx;
        xc = temp.xc;
        YpontRadSx=repmat(YpontRadSx,1,nlay*ps);
        
        RotName   = {'rotor'};
        ShaftName = {'shaft'};
        kkk=1;
        
        % Assign label names
        if (geo.Br==0)
            if isempty(xc)
                BarName=[];
            else
                for kk=1:length(xc)*ps
                    if (YpontRadSx(kk)~=0)
                        num_bar_element=2;
                    else
                        num_bar_element=1;
                    end
                    for jk=1:num_bar_element
                        if (kk<=nlay)
                            ii=kk;
                            ij=jk;
                        else
                            ii=kk-nlay;
                            ij=num_bar_element+jk;
                        end
                        BarName{kkk}={['Air_Bar_',num2str(ii),'_',num2str(ij)]};
                        kkk=kkk+1;
                    end
                end
            end
        else
            BarName=[];
            for kk=1:length(xc)*2*ps
                BarName{kk}={['Plasto_Bar_',num2str(kk)]};
            end
        end
end
