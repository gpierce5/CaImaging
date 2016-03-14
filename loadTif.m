function [ tif ] = loadTif( fpath )
%loadTif imports a tif stack
%   Find and load the motion-corrected tiff stack from 2p
%   Adapted from Clay's impTifStack

%Select folder of ImageJ results files for this mouse
% start_path='C:\Users\Georgia\Data';
% fpath = [uigetdir(start_path) '\'];

%load the tiff stack
filelist=dir([fpath '*_corr.tif']);
tifpath=[fpath filelist.name];
info=imfinfo(tifpath);
nframes=length(info);
%frameWidth = info(1).Width; % width of the first frame (and all others)
%frameHeight = info(1).Height;  % height of the first frame (and all others)

clear stackInfo; % clear this because it might be big

for i=1:nframes
    frame = imread(tifpath, 'tif', i);
    tif(:,:,i) = frame;
end


end

