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

% syreToDxf.m - exports a fem model created by syre to dxf
% input: motorname.mat (created by syre along with motorname.fem)
% output: motorname.dxf, into the new folder motorname\DXF 

clear all; close all; clc;

current_path=cd;
[pathstr, name, ext] = fileparts(current_path);
addpath dxf_conv_fun

load ultimo.mat;

[filemot, pathname] = uigetfile([pathname '\*.mat'], 'Pick a motor');
save ultimo.mat pathname -append;
load([pathname filemot]);

fem=dimMesh(geo,'singt');

if (size(geo.avv,2)<geo.Qs)
    avv=[geo.avv,-geo.avv];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recall nodes, lines, arcs
pathname_DXF=[pathname,filemot(1:end-4),'\DXF\'];

[rotor,BLKLABELSrot,geo] = ROTmatr(geo,fem,mat);
[geo,statore,BLKLABELSstat] = STATmatr(geo,fem);

%%%%%%%%
% export to dxf
mkdir(pathname_DXF);
nomefile=[pathname_DXF,filemot(1:end-4),'.dxf'];
raggi=[];
avvolgimento=[];
magneti=[];

DXFconv(raggi,avvolgimento,rotor,statore,magneti,nomefile);
copyfile([pathname,filemot(1:end-4),'.mat'],[pathname_DXF,filemot(1:end-4),'.mat']);



