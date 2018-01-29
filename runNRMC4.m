function [ timetaken ] = runNRMC
%runNRMC Runs epnev's NoRMCorre on habanero
%INPUTS: name of raw file
%G Pierce Sept 2017
fpath = '/rigel/zi/users/gmp2139';
cd(fpath)
disp('running script')
p = genpath('~/NoRMCorre-master');
addpath(p);
disp('added path')
% fp = '/rigel/users/gmp2139';
% fpath = fullfile(fp,fname);

%filelist = dir(fullfile(fpath,'*.raw'));
%disp('got file')
%for f= 1:length(filelist)
    disp('running through files')
    %fname=filelist(f).name;
    fname = 'rawimage2.raw'
    options_nonrigid = NoRMCorreSetParms('d1',512,'d2',512',...
        'grid_size',[128, 128],'mot_uf',8,'bin_width',50,'max_shift',25,...
        'max_dev',15,'us_fac',50,'phase_flag',true,'iter',1,...
        'output_type','memmap','use_parallel',true,'mem_filename','motion_corrected2.mat');
    disp('made parameters')
    parpool('local',16)
    %disp('got pool')
    %gcp
    tic;
    normcorre_batch(fname,options_nonrigid);
    timetaken = toc
    disp('ran normcorre')
    save('timetaken.mat','timetaken')
    disp('saved')
    
    %NR = matfile('motion_corrected.mat');

%select random set of frames to check by eye (just going to say all vid
%have more than 45000 frames for now
% randFrame = randi(44000);
% saveastiff(uint16(single(NR.Y(:,:,randFrame:(randFrame+999)))),[fullfile(fpath,['vid-frame' num2str(randFrame)]) '.tif']);
    
%end


end

