
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

Content of the folder syreManipulateMM:
scripts for manipulation of the flux linkage tables of synchronous machines

Suggested order of execution:
C_MMLut.m -> C_MtpaMtpvLut.m -> C_OperatingLimits.m -> MaxTw.m

C_MMLut.m
Plots the flux linkage luts in various forms. Prints a text file that can be useful for control. 

input:    1) fdfq_idiq_n256 
          2) ReadParameters.m (see mot_07 example in this folder)
Output:
		  1) figures - Plots the flux linkage luts in various forms.
		  2) FluxTables.txt - Fd Fq tables in the form of txt printed LUTs

C_MtpaMtpvLut.m
Plots the MTPA and MTPV (if applies) curves and prints a text file with MTPA tables (for real time control implementation).

input:    1) fdfq_idiq_n256 
          2) ReadParameters.m
output:   1) AOA\ktMax_idiq.mat (MTPA)
		  2) AOA\kvMax_idiq.mat (MTPV)
          3) tables printed in MtpaTables.txt

C_OperatingLimits.m
Torque and power versus speed profiles, given converter limits Imax and Vmax

input:
          1) fdfq_idiq_n256 
          2) ReadParameters.m
          3) ktMax_idiq.mat and kvMax_idiq.mat

output:
		  1) Plim.mat
		  2) figures saved into a new subfolder

