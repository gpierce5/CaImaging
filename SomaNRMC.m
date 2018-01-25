
%read files from data table
datpath = '/mnt/nas2/homes/georgia/Data/';
%dat = readtable(fullfile(datpath,'DataInfo-S1Imaging1710.xlsx'));

%D = dat(logical(dat.TBMC),:);

%temp just run this one file
fpath = '/mnt/nas2/homes/georgia/Data/SC01-1/171027-001';
fname='Substack2500-5000.tif';
mouse='SC01-1';



for f=1:1%sum(dat.TBMC)
    
    % get the files for this session
    
    %make a big loop
    %fpath = fullfile('/mnt/nas2/homes/georgia/Data/',D.mouse{f},...
        %[num2str(D.day(f)) '-00' num2str(D.nfile(f))]);
    cd(fpath)
    %mouse = D.mouse{f};
    %fname = [D.tiffstack{f} '.tif'];
    tifpath=fullfile(fpath,fname);
    memfilename ='motion_corrected_short.h5';
    
    %set frame to save test video
    %randFrame = randi(44000);
    %numFrames = 1999;
    randFrame=20;
    numFrames=200;
    
    %sesh = strtok(fname,'.');

try
    % set parameters (first try out rigid motion correction)
    
    options_nonrigid = NoRMCorreSetParms('d1',512,'d2',512,...
        'init_batch',200,...
        'bin_width',100,...
        'buffer_width',70,...
        'grid_size',[128, 128],...
        'overlap_pre',[32,32],...
        'mot_uf',8,...
        'max_shift',[50,50],...
        'max_dev',[6,6],...
        'us_fac',50, ...
        'use_parallel', true,...
        'output_type','hdf5','h5_filename',memfilename);
    
    % perform motion correction
    gcp
    tic; [MC,shifts1,template1] = normcorre_batch(tifpath,options_nonrigid); toc
    
    % save tiff of corrected substack
    %NR = matfile(memfilename);
%     tifn=strtok(fname,'.');
%     saveastiff(uint16(single(MC.Y(:,:,randFrame:numFrames+randFrame))),fullfile(fpath,[ tifn '-MC-' num2str(randFrame) '.tif']));
%     
    % catch any errors and text
    load('/home/georgia/phone.mat')
    textme(phone,['Motion correction for ' mouse ' ' fname ' complete!'])
    
catch ME
    load('/home/georgia/phone.mat')
    textme(phone,[':( Error ' ME.identifier ])
    rethrow(ME)
end

end %end file list loop