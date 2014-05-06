% add the required directories to the path
thisfilepath = which('data0.m');

if isoctave
    thisfilepath = fileparts(canonicalize_file_name(thisfilepath));
end
addpath (fullfile (thisfilepath));
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