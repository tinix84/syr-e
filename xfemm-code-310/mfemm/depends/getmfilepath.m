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

function filepath = getmfilepath(mfile)
% getmfilepath: gets the directory containing an mfile
%
% Syntax
%
% filepath = getmfilepath(mfile)
%
% Input
%
% mfile is a string containing the name of the mfile for which the location
% is to be found, the .m extension is optional
% 
% Output
%
% filepath is the directory path containing the specified mfile
%

    if exist(mfile, 'file') == 2
        
       filepath = fileparts(which(mfile));
       
    else
        error('UTILS:nofile', 'm-file does not appear to exist')
    end
    
end