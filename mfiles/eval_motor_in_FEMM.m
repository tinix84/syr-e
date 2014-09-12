% Copyright 2014
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

%% eval_motor_in_FEMM.m
% runs a set of FEA simulations according to the id, iq values given from outside
% - MO_OA is when the script is called by the optimization alogorithm
% - singt and singm are when it is called by

function [out] = eval_motor_in_FEMM(geo,io,gamma_in,eval_type)

switch eval_type
    case 'MO_OA'
        gamma = gamma_in;
        nsim = geo.nsim_MOOA;
        delta_sim = geo.delta_sim_MOOA;
    case 'MO_GA'
        gamma = gamma_in;
        nsim = geo.nsim_MOOA;
        delta_sim = geo.delta_sim_MOOA;
    case 'singt'
        nsim = geo.nsim_singt; delta_sim = geo.delta_sim_singt;gamma = gamma_in;        
end


SOL = [];
T = zeros(1,length(gamma)); fd = T; fq = T; ripple = T;
Pfe_S_tot = T; Pfe_R_tot = T;


kk = 1;
SOL(:,:,kk) = simulate_xdeg(geo,nsim,delta_sim,io,gamma(kk),eval_type); % risultati
SOL(:,2,kk)=SOL(:,2,kk)/geo.Nbob;
SOL(:,3,kk)=SOL(:,3,kk)/geo.Nbob;
SOL(:,4,kk)=SOL(:,4,kk)*geo.Nbob;
SOL(:,5,kk)=SOL(:,5,kk)*geo.Nbob;

ris_sim = SOL(:,:,kk);

T(kk) = abs(mean(ris_sim(:,6)));
ripple(kk) = std(ris_sim(:,6));
fd(kk) = mean(ris_sim(:,4));
fq(kk) = mean(ris_sim(:,5));

out.SOL = SOL(:,:,:);
out.Tn = T;
out.ripple_pu = ripple./T;
out.fd = fd;
out.fq = fq;

save out

