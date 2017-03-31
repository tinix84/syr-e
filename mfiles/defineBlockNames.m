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

function BarName = defineBlockNames(temp,geo,mat)

% define names of rotor barriers (not used in FEMM, needed for
% compatibility with other finite-element codes)

ps = geo.ps;
nlay = geo.nlay;
RotType = geo.RotType;
seg = geo.dx;

if strcmp(RotType,'SPM')
    RotName   = {'rotor'};
    ShaftName = {'shaft'};
    PMName    = {'permanent magnet'};
    kkk=1;
    BarName = [];
    for kk = 1:ps*seg
        BarName{kkk} = {['PM_',num2str(kk)]};
        kkk=kkk+1;
    end
else
    YpontRadSx = temp.YpontRadSx;
    xmag = temp.xmag;
    YpontRadSx=repmat(YpontRadSx,1,nlay*ps);
    
    RotName   = {'rotor'};
    ShaftName = {'shaft'};
    kkk=1;
    
    % Assigment of label name
    if (mat.LayerMag.Br==0)
        if isempty(xmag)
            BarName=[];
        else
            for kk=1:nlay*ps
                if (YpontRadSx(kk)~=0)
                    num_bar_element=4;
                else
                    num_bar_element=2;
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
        for kk=1:length(xmag)*2*ps
            BarName{kk}={['Plasto_Bar_',num2str(kk)]};
        end
    end
end
