foldername = '/mnt/nas2/homes/georgia/Data/SC01-1/171027-001/';
        % folder where all the files are located. Currently supported .tif,
        % .hdf5, .raw, .avi, and .mat files
%files = subdir(fullfile(foldername,'*.tif'));   % list of filenames (will search all subdirectories)
FOV = [512 512];%size(read_file(files(1).name,1,1));
numFiles = 1; %length(files);
fpath = fullfile(foldername,'motion_corrected.h5');
motion_correct = true;

%% downsample mat files and save into a single memory mapped matlab file
disp('downsampling')
if motion_correct
    h5_files = subdir(fullfile(foldername,'motion_corrected.h5'));  % list of h5 files (modify this to list all the motion corrected files you need to process)
else
    h5_files = subdir(fullfile(foldername,'*_mc.h5'));
end

fr = 30;                                         % frame rate
tsub = 7;                                        % degree of downsampling (for 30Hz imaging rate you can try also larger, e.g. 8-10)
ds_filename = fullfile(foldername,'ds_data.mat');
data_type = class(read_file(fpath,1,1));
data = matfile(ds_filename,'Writable',true);
data.Y  = zeros([FOV,0],data_type);
data.Yr = zeros([prod(FOV),0],data_type);
data.sizY = [FOV,0];
F_dark = Inf;                                    % dark fluorescence (min of all data)
batch_size = 1000;                          % read chunks of that size
batch_size = round(batch_size/tsub)*tsub;        % make sure batch_size is divisble by tsub
Ts = zeros(numFiles,1);                          % store length of each file
cnt = 0;                                         % number of frames processed so far
tt1 = tic;
for i = 1%:numFiles
    name = h5_files(i).name;
    info = h5info(name);
    dims = info.Datasets.Dataspace.Size;
    ndimsY = length(dims);                       % number of dimensions (data array might be already reshaped)
    Ts(i) = dims(end);
    Ysub = zeros(FOV(1),FOV(2),floor(Ts(i)/tsub),data_type);
    data.Y(FOV(1),FOV(2),sum(floor(Ts/tsub))) = zeros(1,data_type);
    data.Yr(prod(FOV),sum(floor(Ts/tsub))) = zeros(1,data_type);
    cnt_sub = 0;
    for t = 1:batch_size:Ts(i)
        Y = bigread2(name,t,min(batch_size,Ts(i)-t+1));    
        F_dark = min(nanmin(Y(:)),F_dark);
        ln = size(Y,ndimsY);
        Y = reshape(Y,[FOV,ln]);
        Y = cast(downsample_data(Y,'time',tsub),data_type);
        ln = size(Y,3);
        Ysub(:,:,cnt_sub+1:cnt_sub+ln) = Y;
        cnt_sub = cnt_sub + ln;
    end
    data.Y(:,:,cnt+1:cnt+cnt_sub) = Ysub;
    data.Yr(:,cnt+1:cnt+cnt_sub) = reshape(Ysub,[],cnt_sub);
    toc(tt1);
    cnt = cnt + cnt_sub;
    data.sizY(1,3) = cnt;
end
data.F_dark = F_dark;

%% now run CNMF on patches on the downsampled file, set parameters first
disp('setting params')
sizY = data.sizY;                       % size of data matrix
patch_size = [128,128];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [32,32];                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 250;                                            % number of components to be found
tau = 7;                                          % std of gaussian kernel (half size of neuron) 
p = 2;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.7;                                  % merging threshold
sizY = data.sizY;


options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'deconv_method','constrained_foopsi',...    % neural activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'ssub',10,...                                % spatial downsampling when processing
    'tsub',4,...                                % further temporal downsampling when processing
    'merge_thr',merge_thr,...                   % merging threshold
    'gSig',tau,... 
    'max_size_thr',400,'min_size_thr',40,...    % max/min acceptable size for each component
    'spatial_method','regularized',...          % method for updating spatial components
    'se', strel('disk',3,0),... % morphological element for dilation (default: strel('disk',1,0))
    'df_prctile',50,...                         % take the median of background fluorescence to compute baseline fluorescence 
    'fr',fr/tsub,...                            % downsamples
    'space_thresh',0.5,...                      % space correlation acceptance threshold
    'min_SNR',2.0,...                           % trace SNR acceptance threshold
    'cnn_thr',0.2,...                           % cnn classifier acceptance threshold
    'nb',1,...                                  % number of background components per patch
    'gnb',3,...                                 % number of global background components
    'decay_time',0.5...                        % length of typical transient for the indicator used
    );

%% Run on patches (the main work is done here)
disp('running CNMF on patches')

[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,0,options);  % do not perform deconvolution here since
                                                          % we are operating on downsampled data         
                                                            
%% compute correlation image on a small sample of the data (optional - for visualization purposes) 
disp('computing correlation image')
Cn = correlation_image_max(data.Y,8);

%% classify components
disp('classifying components')
rval_space = classify_comp_corr(data,A,C,b,f,options);
ind_corr = rval_space > options.space_thresh;           % components that pass the correlation test
                                        % this test will keep processes
                                        
                                        
                                        
%% further classification with cnn_classifier
disp('even more classification')
try  % matlab 2017b or later is needed
    [ind_cnn,value] = cnn_classifier(A,FOV,'/home/georgia/Documents/MATLAB/cnn_model.h5',options.cnn_thr);
catch
    ind_cnn = true(size(A,2),1);                        % components that pass the CNN classifier
end     
                            
%% event exceptionality
disp('event exceptionality??')
fitness = compute_event_exceptionality(C+YrA,options.N_samples_exc,options.robust_std);
ind_exc = (fitness < options.min_fitness);

%% select components
disp('selecting components')
keep = (ind_corr | ind_cnn) & ind_exc;

%% run GUI for modifying component selection (optional, close twice to save values)
% run_GUI = false; %can't get the GUI to work
% if run_GUI
%     Coor = plot_contours(A,Cn,options,1); close;
%     %GUIout = ROI_GUI(Y,A,P,options,Cn,C,b,f); %original demo_GUI statement
%     %GUIout = ROI_GUI(Ysub,A,P,options,Cn,C,b,f); %gmp can open it with
%     %this at least
%     GUIout = ROI_GUI(A,options,Cn,Coor,keep,ROIvars); %original
%     %run_pipeline statment
%     options2 = GUIout{2};
%     keep2 = GUIout{3};    
% end

%% view contour plots of selected and rejected components (optional)
disp('viewing contour plots, showing active components as keep')
throw = ~keep;
Coor_k = [];
Coor_t = [];
figure; 
    %without numbers
    ax1 = subplot(121); plot_contours(A(:,keep),Cn,options,0,[],Coor_k,[],1,find(keep)); title('Selected components','fontweight','bold','fontsize',14);
    ax2 = subplot(122); plot_contours(A(:,throw),Cn,options,0,[],Coor_t,[],1,find(throw));title('Rejected components','fontweight','bold','fontsize',14);
    linkaxes([ax1,ax2],'xy')
    colormap gray
    
%% keep only the active components    
  
A_keep = A(:,keep);
C_keep = C(keep,:);

%% manually accept/reject components
 %to view component numbers if that line says i=indshow not i=1:ind_show
    %ax1 = subplot(121); plot_contours(A(:,keep),Cn,options,1,[],Coor_k,[],length(find(keep))); title('Selected components','fontweight','bold','fontsize',14);
figure
%inputs:
%plot_contours(Aor,Cn,options,display_numbers,max_number,Coor, ln_cl, ind_show,cm)
compind = find(keep);
manualkeep = keep;
for curr = 1:sum(keep)
plot_contours(A_keep,Cn,options,1,curr,Coor_k,[],curr,[]); title(['component ' num2str(curr) ' of ' num2str(sum(keep))],'fontweight','bold','fontsize',14);
colormap gray
response = input('Keep (1) or Reject (0) component?');
if response==1 || response==0
manualkeep(compind(curr)) = response;
else
    disp('you did not answer 0 or 1, try again next time.')
end

end


%% deconvolve (downsampled) temporal components plot GUI with components (optional)

% tic;
% [C_keep,f_keep,Pk,Sk,YrAk] = update_temporal_components_fast(data,A_keep,b,C_keep,f,P,options);
% toc
% 
% plot_components_GUI(data,A_keep,C_keep,b,f,Cn,options)
 if exist('YrAk','var'); R_keep = YrAk; else; R_keep = YrA(keep,:); end
    
%% extract fluorescence on native temporal resolution

options.fr = options.fr*tsub;                   % revert to origingal frame rate
N = size(C_keep,1);                             % total number of components
T = sum(Ts);                                    % total number of timesteps
C_full = imresize(C_keep,[N,T]);                % upsample to original frame rate
R_full = imresize(R_keep,[N,T]);                % upsample to original frame rate
F_full = C_full + R_full;                       % full fluorescence
f_full = imresize(f,[size(f,1),T]);             % upsample temporal background

S_full = zeros(N,T);

P.p = 0;
ind_T = [0;cumsum(Ts(:))];
options.nb = options.gnb;
for i = 1:numFiles
    inds = ind_T(i)+1:ind_T(i+1);   % indeces of file i to be updated
    [C_full(:,inds),f_full(:,inds),~,~,R_full(:,inds)] = update_temporal_components_fast(h5_files(i).name,A_keep,b,C_full(:,inds),f_full(:,inds),P,options);
    disp(['Extracting raw fluorescence at native frame rate. File ',num2str(i),' out of ',num2str(numFiles),' finished processing.'])
end

%% extract DF/F

[F_dff,F0] = detrend_df_f(A_keep,[b,ones(prod(FOV),1)],C_full,[f_full;-double(F_dark)*ones(1,T)],R_full,options);
%F_dff is detrended fluor in DF/F, F0 is baseline fluor for each trace