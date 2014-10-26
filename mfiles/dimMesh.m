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

% adapts
function fem = dimMesh(geo,eval_type)

res_max=1/3*geo.g;
% La dimensione degli elementi della mesh al traferro deve essere un
% multiplo (o sottomultiplo) del passo di simulazione.
if strcmp(eval_type,'MO_GA')||strcmp(eval_type,'MO_OA')
    fem.res_traf=1/geo.p;
    fem.res=geo.K_mesh_MOOA*fem.res_traf;
else
    %     res_multiplo=(geo.delta_sim_singt/(geo.nsim_singt)*180/pi)*...
    %         (geo.r+0.5*geo.g);
    %     if res_multiplo<1.5*res_max
    %         fem.res_traf=res_multiplo;
    %     else
    %         % Se res_multiplo e' maggiore del massimo (res_max), impongo che il
    %         % valore di risoluzione scelto per il traferro sia un sottomultiplo di
    %         % res_multiplo
    %         cerco_res=res_multiplo;
    %         divido=2;
    %         while cerco_res>1.5*res_max
    %             cerco_res=res_multiplo/divido;
    %             divido=divido+1;
    %         end
    %         fem.res_traf=cerco_res;
    %     end
    %     % - Resto della macchina
    %     fem.res=1.5*geo.K_mesh*fem.res_traf;
    
    fem.res_traf=1/geo.p * geo.K_mesh/geo.K_mesh_MOOA;
    fem.res=geo.K_mesh_MOOA*fem.res_traf;
    
end

