%% read big tiff
fpath = '/home/georgia/ToBeMotionCorrected/9029-3';
fname = '170412-4x_00001.tif';
mouse = '9029-3';
sesh = strtok(fname,'.');

cd(fpath)
Y = ScanImageTiffReader(fname).data();
Y = double(Y);

%% compute metrics for original
nnY = quantile(Y(:),0.005);
mmY = quantile(Y(:),0.995);

%consider adding portion for increasing template init batch if the first
%section has a weird correlation

[cY,mY,vY] = motion_metrics(Y,10,10);
d = ndims(Y);
T = length(cY);

%% perform non-rigid motion correction
options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[128,128],...
    'mot_uf',8,'bin_width',50,'max_shift',25,'max_dev',10,'us_fac',50,...
    'phase_flag',true,'iter',1,...
    'memmap',true,'mem_filename',[sesh '_MC.mat'],'use_parallel',true);

tic; 
[NR,shifts1,template1] = normcorre_batch(Y,options_nonrigid); 
toc

%% compute metrics for NR
[cNR,mNR,vNR] = motion_metrics(NR,10,10);

%% plot metrics
figure;
    ax1 = subplot(2,2,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off;
        title('mean raw data','fontsize',14,'fontweight','bold')
    ax2 = subplot(2,2,2); imagesc(mNR,[nnY,mmY]);  axis equal; axis tight; axis off;
        title('mean raw data','fontsize',14,'fontweight','bold')
    linkaxes([ax1,ax2],'xy')
    
    subplot(2,2,3); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); 
        title('correlation coefficients','fontsize',14,'fontweight','bold')
    subplot(2,2,4); scatter(cY,cNR); hold on; plot([0.9*min(cY),1.05*max(cNR)],[0.9*min(cY),1.05*max(cNR)],'--r'); axis square;
        xlabel('raw data','fontsize',14,'fontweight','bold'); 
        ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');
    
%% plot shifts        
% shifts_r = horzcat(shifts1(:).shifts)';
% shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
% shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
% shifts_x = squeeze(shifts_nr(:,1,:))';
% shifts_y = squeeze(shifts_nr(:,2,:))';
% 
% patch_id = 1:size(shifts_x,2);
% str = strtrim(cellstr(int2str(patch_id.')));
% str = cellfun(@(x) ['patch # ',x],str,'un',0);
% 
% figure;
%     ax1 = subplot(311); plot(1:T,cY,1:T,cM2); legend('raw data','non-rigid'); 
%         title('correlation coefficients','fontsize',14,'fontweight','bold')
%             set(gca,'Xtick',[])
%     ax2 = subplot(312); plot(shifts_x); hold on; plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
%             set(gca,'Xtick',[])
%     ax3 = subplot(313); plot(shifts_y); hold on; plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
%             xlabel('timestep','fontsize',14,'fontweight','bold')
%     linkaxes([ax1,ax2,ax3],'x')

%% save a simple tiff for checking in imagej

savemehere = fullfile(fpath,sesh);
% tifoptions.color = false;
% tifoptions.compress = 'no';
% tifoptions.message = true;
% tifoptions.append = false;
% tifoptions.overwrite = false;
% tifoptions.big = false;

saveastiff(uint16(single(Y(:,:,1:1000))),[savemehere '_MCtest-beginRAW.tif'])%,tifoptions)

if mod(T,10)==0
    ct=1;
    f1=1;
    while f1<size(Y,3)
        f2=f1+9;
        Z(:,:,ct) = nanmean(NR(:,:,f1:f2),3);
        ct=ct+1;
        f1=f2+1;
    end
    saveastiff(uint16(single(Z)),[savemehere '_MCtest.tif'])
    
    ct=1;
    f1=1;
    while f1<size(Y,3)
        f2=f1+9;
        Z_raw(:,:,ct) = mean(Y(:,:,f1:f2),3);
        ct=ct+1;
        f1=f2+1;
    end
    saveastiff(uint16(single(Z_raw)),[savemehere '_MCtest_RAW.tif'])
else
    disp('mod issue')
end

