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

%% MG 2013/07/02 
%% Export MOGA to DXF 
% Il seguente scrip realizza l'export tra la geometria di MOGA e DXF

%%%%%%%%%%%%%%%
ROTmatr;
%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
STATmatr
%%%%%%%%%%%%%%%%
% closefemm;

%%  %%%%%
%% Salvataggio del file dxf
%%  %%%%%
mkdir(pathname_DXF);
nomefile=[pathname_DXF,filemot(1:end-4),'.dxf'];

raggi=[];
avvolgimento=[];
magneti=[];
rotore=[];

%% Creazione e salvataggio del file dxf:
DXFconv(raggi,avvolgimento,rotore2,statore,magneti,nomefile);







