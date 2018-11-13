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

function BarName = defineBlockNames(BarCenter,geo)

% define names of rotor barriers (not used in FEMM, needed for
% compatibility with other finite-element codes)

%% Main revision
% 01/08/2018 (Simone Ferrari): now BarName is write starting from
% BarCenter. This method is more robust respect the previous one and
% doesn't need code changes if the geometry is changed.

indAir=1;
indPMs=1;
ind=1;

rotType=geo.RotType;

% Material codes: refer to BLKLABELS.materials positions
codMatFe    = 5;
codMatBar   = 6;
codMatShaft = 7;
codMatPM    = 6;
codMatAir   = 1;

for ii=1:length(BarCenter(:,3))
    if (BarCenter(ii,3)==codMatAir && ~strcmp(rotType,'SPM'))
        BarName{ind}=['Air_Bar_' int2str(indAir)];
        ind=ind+1;
        indAir=indAir+1;
    elseif BarCenter(ii,3)==codMatPM
        BarName{ind}=['Magnet_Bar_' int2str(indPMs)];
        ind=ind+1;
        indPMs=indPMs+1;
    end
end



%% previous code, changed on 01/08/2018 by Simone Ferrari
% function BarName = defineBlockNames(temp,geo,mat)
% 
% ps = geo.ps;
% nlay = geo.nlay;
% RotType = geo.RotType;
% seg = geo.dx;
% 
% if strcmp(RotType,'SPM')
%     RotName   = {'rotor'};
%     ShaftName = {'shaft'};
%     PMName    = {'permanent magnet'};
%     kkk=1;
%     BarName = [];
%     for kk = 1:ps*seg
%         BarName{kkk} = {['PM_',num2str(kk)]};
%         kkk=kkk+1;
%     end
% else
%     YpontRadSx = temp.YpontRadSx;
%     xmag = temp.xmag;
%     YpontRadSx=repmat(YpontRadSx,1,nlay*ps);
%     
%     RotName   = {'rotor'};
%     ShaftName = {'shaft'};
%     kkk=1;
%     
%     % Assigment of label name
%     if (mat.LayerMag.Br==0)
%         if isempty(xmag)
%             BarName=[];
%         else
%             for kk=1:nlay*ps
%                 if (YpontRadSx(kk)~=0)
%                     num_bar_element=4;
%                 else
%                     num_bar_element=2;
%                 end
%                 for jk=1:num_bar_element
%                     if (kk<=nlay)
%                         ii=kk;
%                         ij=jk;
%                     else
%                         ii=kk-nlay;
%                         ij=num_bar_element+jk;
%                     end
%                     BarName{kkk}={['Air_Bar_',num2str(ii),'_',num2str(ij)]};
%                     kkk=kkk+1;
%                 end
%             end
%         end
%     else
%         BarName=[];
%         for kk=1:length(xmag)*2*ps
%             BarName{kk}={['Plasto_Bar_',num2str(kk)]};
%         end
%     end    
% end
