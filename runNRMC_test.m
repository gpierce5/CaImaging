function [ ] = runNRMC
%runNRMC Runs epnev's NoRMCorre on habanero
%INPUTS: name of raw file
%G Pierce Sept 2017
fpath = '/rigel/zi/users/gmp2139';
cd(fpath)
p = genpath('~/NoRMCorre-master');
addpath(p);
    fname = 'rawimage.raw';
    options_nonrigid = NoRMCorreSetParms('d1',512,'d2',512',...
        'grid_size',[128, 128],'mot_uf',8,'bin_width',50,'max_shift',25,...
        'max_dev',15,'us_fac',50,'phase_flag',true,'iter',1,...
        'output_type','memmap','use_parallel',true);
    disp('made parameters')
    parpool('local',24)
    normcorre_batch(fname,options_nonrigid);
end

