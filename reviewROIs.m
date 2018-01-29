function [ROIkeep] = reviewROIs(cnmfpath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%   INPUT:  cnmfpath - path to file of CNMF output



load(cnmfpath)
%works only for windows looking at nas2...
[~,remnant1] = strtok(cnmfpath,'\');
[~,remnant2] = strtok(remnant1,'\');
[~,remnant3] = strtok(remnant2,'\');
[~,remnant4] = strtok(remnant3,'\');
[mouse,remnant5]=strtok(remnant4,'\');
[sesh,~] = strtok(remnant5,'\');
disp(['Reviewing ' mouse ' ' sesh])

% keep = (ind_corr | ind_cnn) & ind_exc;
% A_keep = A(:,keep);
% C_keep = C(keep,:);

if ~exist('F_dff','var')
    FOV = [512 512]
    [F_dff,F0] = detrend_df_f(A_keep,[b,ones(prod(FOV),1)],C_full,[f_full;-double(F_dark)*ones(1,T)],R_full,options);
end

cafile = [sesh '_CaTraces.mat'];
if exist(cafile,'file')==2 %amanda's nam2
    load(cafile)
else
    ROIkeep = zeros(1,size(A_keep,2));
    save(cafile, 'ROIkeep')
end

if ~exist('sframe','var'), sframe = 1;end
if ~exist('num2read','var'),num2read = 9000;end

%T = size(F_dff,2); already saved T
d1 = dims(1); %sqrt(size(P.sn,1));%size(Yr,1);
d2 = dims(2); %d1
%Cn =  reshape(P.sn,d1,d2); %correlation_image(Y); %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
C_df = F_dff; %not sure if C_full is supposed to be the new C_df... %full(C_df);
%I think we want the delta f calculation here...

%Display background
j = figure('Position',[100 100 800 800]);
set(j,'WindowStyle','docked');
yMax = quantile(C_df(:),0.9990); if issparse(yMax),yMax = full(yMax);end
yMin = quantile(C_df(:),0.0010); if issparse(yMin),yMin = full(yMin);end
subplot(4,2,1:6)
imagesc(Cn), colormap('gray'),axis square, hold on
subplot(4,2,7:8); xlabel('frame #'); ylabel('df / F'), axis([0 num2read yMin yMax])
plot(1:num2read,C_df(:,1:num2read),'Color',[0.5 0.5 0.5]),hold on

%Check if have previously saved some ROIs. If so, plot in red/green
if isequal(nansum(ROIkeep),0)
    start = 1;
else
    b = input('start at beginning of ROIs? ');
    if isequal(b,'y')
        start = 1;
        ROIkeep = zeros(1,size(A_keep,2)-1);
    else
        start = find(ROIkeep == 1,1,'last');
        for i = 1:start
            figure(j)
            subplot(4,2,1:6)
            A_temp = full(reshape(A_keep(:,i),d1,d2));
            if ROIkeep(i) == 1; plotSingleROI(A_temp,options.maxthr,'g');
            else plotSingleROI(A_temp,options.maxthr,'r'); end
        end
    end
end

%amanda eliminates non-active ROIs like this
%probably could use exceptionality indices from epnev here
% m = nanmean(C_df,2);
% s = max(C_df')';
% s2 = std(C_df')';
% ROIkeep(s-s2 <= 1.4) = NaN; %Automatically discard any flat traces
% ROIkeep(end) = [];
%
% A2(:,isnan(ROIkeep)) = NaN;
% A_or(:,isnan(ROIkeep)) = NaN;
% C2(isnan(ROIkeep),:) = NaN;
% C_df(isnan(ROIkeep),:) = NaN;

for i = start:size(A_keep,2) %loop through each ROI
    if ~isnan(ROIkeep(i))
        figure(j)
        subplot(4,2,1:6)
        if i > 1
            %figure(j)
            %subplot(2,2,1)
            %imagesc(Cn), colormap('gray'),axis square
            %for jj = 1:i-1
            jj=i-1;
                A_temp = full(reshape(A_keep(:,jj),d1,d2));
                if ROIkeep(jj) == 1; plotSingleROI(A_temp,options.maxthr,jj,'g'); else plotSingleROI(A_temp,options.maxthr,'r'); end
            %end
        end
        A_temp = full(reshape(A_keep(:,i),d1,d2));
        roi = plotSingleROI(A_temp,options.maxthr,'b');
      
%         subplot(2,2,2)
%         cla
%         dispROI(A_or,d1,d2,i,T,options);
        
        subplot(4,2,7:8); xlabel('frame #'); ylabel('df / F'), axis([0 num2read yMin yMax])
        %cla
%         plot(1:size(C_df,2)-5,C_df(:,1:end-5),'Color',[0.5 0.5 0.5]),hold on
%         plot(1:size(C_df,2)-5,C_df(i,1:end-5),'b'), hold off
       roidff= plot(1:num2read,C_df(i,1:num2read),'b');
        
        %axis tight
        drawnow limitrate
        
        %figure(j)
        r = [];
        if ~isnan(ROIkeep(i))
            while isempty(r)
                r = input('Keep ROI? (1 = yes, 0 = no): ');
                if ~abs(r) > 1
                    r = [];
                end
                
            end
            if r == 0 || r==1
            ROIkeep(i) = r;
            delete(roi);
            delete(roidff);
            elseif r == 'q'
                break;
            end
        else
            disp('Automatically rejecting ROI...')
            ROIkeep(i) = 0;
                        delete(roi);
            delete(roidff);
            %pause(0.5)
        end
%         if isequal(r,'n')
%             break;
%         end
        save(cafile,'ROIkeep','C_df','-append')
    end
end










end

