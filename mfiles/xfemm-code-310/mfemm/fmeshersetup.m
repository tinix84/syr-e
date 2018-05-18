function fmeshersetup(dodebug, verbose)
% compiles the fmesher mexfunction 

% Copyright 2012 Richard Crozier
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

    if nargin < 1
        dodebug = true;
    end
    
    if nargin < 2
        verbose = true;
    end    

    % store the current directory
    origdir = pwd;
    
    % return to original dir on interruption or completion
    %OC = onCleanup (@() cd(origdir));
    
    % change to the mfemm directory (the directory this file is in)
    cd(fileparts(which('fmeshersetup.m')));
    
    % make architecture specific mex directory if it doesn't already exist
    warning off MATLAB:MKDIR:DirectoryExists
    mexdir = ['xfemm_mex_files_for_' computer('arch')];
    mkdir (mexdir);
    warning on MATLAB:MKDIR:DirectoryExists

    cd (mexdir);
    
    % set some common compiler flags, we replace all calls to printf to
    % calls to mexPrintf
    if dodebug
        common_compiler_flags = {'-g'};
    else
        common_compiler_flags = {}; %common_compiler_flags = {'-D"PRINTF=mexPrintf"'};
    end
   
    libcommands = {'-I"../pfemm/fmesher"', ...
                   '-I"../pfemm/libfemm"', ...
                   '-I"../pfemm/libfemm/liblua"', ...
                   '../pfemm/fmesher/libfmesher.a'};

    % add an appropriate flag for the trilibrary C external depending on
    % computer architecture
    if ispc
        trilibraryflag = {'-DCPU86'};
    else
        trilibraryflag = {'-DLINUX'};
    end
    
    if verbose
        common_compiler_flags = [common_compiler_flags, {'-v'}];
    end

    % construnct the command string 
    mexcommands = [ common_compiler_flags, ...
                    trilibraryflag, ...
                    { '../mexfmesher.cpp' } ...
                    libcommands ...
                  ];

    % call mex with the appropriately constructed commands
    if isoctave
        mkoctfile('--mex', mexcommands{:});
    else
        mex(mexcommands{:});
    end
    cd(origdir)
end
