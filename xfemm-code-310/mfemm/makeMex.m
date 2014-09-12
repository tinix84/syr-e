mexdir = fullfile(['xfemm_mex_files_for_' computer('arch')]);
eval(['cd ' mexdir])
delete *mex*
eval(['cd ..'])
fmeshersetup
fsolversetup
fpprocsetup
