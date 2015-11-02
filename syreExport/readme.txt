
% Copyright 2015
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Content of the folder syreExport:
scripts for export syre to dxf and vice versa, scripts for export dxf to Magnet and run Magnet automatically

syreToDxf.m - exports a fem model created by syre to dxf
% input: motorname.mat (created by syre along with motorname.fem)
% output: motorname.dxf, into the new folder motorname\DXF 

BuildMachine.m - builds the machine model into Magnet by Infolytica
% input:  motorname.mat (made by syre) and motorname.dxf (from syreToDxf.m)
% output: motorname.mn in the same folder

RUN_SIM.m - Simulates single id,iq,rpm operating points or the entire flux linkage map (similar to post processing in syre).

% input: motorname.mat (made by syre) and motorname.mn (from BuildMachine.m)
% output: mat files in new subfolder CaseBackup

% Select SIMULATION or IDENTIFICATION modes, accordingly.
% SIMULATION stands for a single or multiple id, iq combination, IDENTIFICATION is same as gamma = 1000 in syre post proc. 
% You can use RUN_SIM.m to replicate the magnetic model done with FEMM (flux curves must be the same) and, more importantly, to evaluate iron and PM loss of your machines.

% ELAB_SIM.m - Processes data from single or multiple SIMULATION from RUN_SIM