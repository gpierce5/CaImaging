function [ roi ] = plotSingleROI( A_temp,thr, ln_cl )
%plotSingleROI Plot specified ROI i. Does not plot background image
%   Inputs: Aor - ordered A matrix
%       thr - threshold from options.maxthr
%       i - index of the ROI to plot
% Pulled from epnev's CaImAn-MATLAB plot_contours


                A_temp = medfilt2(A_temp,[3,3]);
                A_temp(A_temp<thr*max(A_temp(:))) = 0;
                BW = bwareafilt(A_temp>0,1);                
                BW2 = bwboundaries(BW);
                if ~isempty(BW2)
                    for ii = 1:length(BW2)
                        BW2{ii} = fliplr(BW2{ii});
                        roi = plot(BW2{ii}(:,1),BW2{ii}(:,2),'Color',ln_cl, 'linewidth', 2);
                    end
%                     CC{1} = BW2{1}';
%                     fp = find(BW);
%                     [ii,jj] = ind2sub([d1,d2],fp);
%                     CR{1,1} = [ii,jj]';
%                     CR{1,2} = A_temp(fp)';
                end

end

