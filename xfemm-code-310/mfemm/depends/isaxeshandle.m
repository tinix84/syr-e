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

function result = isaxeshandle(h)
% isaxeshandle: returns true if an object is a handle to an axes object
% (and not a legend or colorbar) 

    result = zeros(size(h));
    
    handles = ishandle(h);
    
    result(handles) = strcmp( get(h(handles), 'type'), 'axes' ) & ...
                                ~(strcmp( get(h(handles), 'Tag'), 'legend') | ...
                                  strcmp( get(h(handles), 'Tag'), 'Colorbar'));

end