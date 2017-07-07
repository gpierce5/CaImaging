function [ output_args ] = reduceTiff( fpath,fname,onServer )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%add folder to path
if onServer
    p = genpath(~/ScanImageTiffReader/ScanImageTiffReader-v1.2-pre-4-gb31cc23-Linux/share/matlab);
    addpath(p);
end



end

