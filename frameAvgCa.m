function [ frameAvgDF, frameAvg ] = frameAvgCa( fpath )
%frameAvgCa whole-frame avg fluorescence per frame
%   Take 2p tiff stack/video and output deltaF/F for
%   the session. Optional second variable for looking at unprocessed
%   fluorescence.

tif = loadTif(fpath);

% take average value for each frame
frameAvg = mean(mean(tif,1),2);  
frameAvg = squeeze(frameAvg); %make the variable [number of frames x 1]

%deltaF/F = (F-Fo)/Fo, Fo is the baseline
%Approximate a baseline by taking a running mean, may want to re-evaluate
%whether this is the best method
%note: in runmean, give half the window size because it centers the window on
%either side of the current index
meanbin=2000; %2000 frames is approx 4 min at 128ms per frame, so averaging over 8 min?
Fo = runmean(frameAvg, meanbin);
frameAvgDF = (frameAvg-Fo)./Fo*100;

figure
ax1 = subplot(2,1,1);
plot(frameAvg)
ax2 = subplot(2,1,2);
plot(frameAvgDF,'m')
linkaxes([ax1, ax2],'x')

end

