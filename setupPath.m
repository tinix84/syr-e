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

% add the required directories to the path
thisfilepath = fileparts(which('data0'));

%if isoctave
%    thisfilepath = fileparts(canonicalize_file_name(thisfilepath));
%    endfileName = '\';
%end
addpath (fullfile ([thisfilepath]));
addpath (fullfile (thisfilepath,'dxf_conv_fun'));
addpath (fullfile (thisfilepath,'mfiles'));
addpath (fullfile (thisfilepath,'MODE'));
addpath (fullfile (thisfilepath,'results'));
addpath (fullfile (thisfilepath,'ValutaMOTOR'));
xfemmPath = fullfile (thisfilepath, 'xfemm-code-310','mfemm');
addpath (fullfile (xfemmPath));
addpath (fullfile (xfemmPath, 'preproc'));
addpath (fullfile (xfemmPath, 'postproc'));
addpath (fullfile (xfemmPath, 'examples'));
addpath (fullfile (xfemmPath, 'depends'));
addpath (fullfile (xfemmPath, 'visualisation'));
mexdir = fullfile(xfemmPath, ['xfemm_mex_files_for_' computer('arch')]);
addpath (mexdir);
savepath