function [galvoFrameInd] = getResonantFrames(event)
%getResonantFrames Get frame times (in units of ntrode index) for our resonant scanner situation
%%% INPUT VARIABLES: %%%%
% event: raw signal from ntrode for resonant galvo

%%% OUTPUT VARIABLES: %%%%
%resFrameInd: vector of indices for the last scan of each frame

%took these out, do you still want them?:
%galvoStartInd: index of the start of the galvo signal
%galvoEndInd: index of the end of the galvo signal

[pks,locs] = findpeaks(event,'MinPeakDistance',400);
%keep an eye out if 400 continues to be a good min pk distance
%for now we can use the peak of the galvo signal because it's clean
%enough and because Ca will have a slight lag anyway, so the imaging
%comes after the behavior anyway. We'll take time before the peak
%as events going on in that frame of activity

%now to find the start and the end of the imaging
minpk = min(pks);
maxpk = max(pks);
midpk = minpk/maxpk;
framepks = pks(pks>midpk);
frameTimes = locs(pks>midpk);
%plot(frameTimes,framepks,'kp')

%were there multiple imaging times? Only keep the last one
%seems like big positive going derivatives might be the sign of beginning
%imaging
dEvent=diff(event);
[dpks,dlocs]=findpeaks(dEvent);
startScan = dlocs(dpks>2);
if startScan>1
    disp('Oops, ' num2str(length(startScan)) ' galvo runs. Keeping final session')
    frameTimes=frameTimes(frameTimes >= startScan(end));
end

disp(['found ' num2str(length(frameTimes)) ' frames'])
imgtime = (length(frameTimes)/30)/60; %minutes imaged in resonant
disp(['If 30Hz, looks like you imaged ' num2str(imgtime) ' min'])

galvoFrameInd = frameTimes;

end

