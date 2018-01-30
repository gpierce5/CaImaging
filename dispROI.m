

function dispROI(A_or,d1,d2,i,T,options)
cla
nr = size(A_or,2);
cm = com(A_or,d1,d2);
sx = min([options.sx,floor(d1/2),floor(d2/2)]);
int_x = zeros(nr,2*sx);

Atemp = reshape(A_or(:,i),d1,d2);
int_x(i,:) = round(cm(i,1)) + (-(sx-1):sx);
if int_x(i,1)<1
    int_x(i,:) = int_x(i,:) + 1 - int_x(i,1);
end
if int_x(i,end)>d1
    int_x(i,:) = int_x(i,:) - (int_x(i,end)-d1);
end
int_y(i,:) = round(cm(i,2)) + (-(sx-1):sx);
if int_y(i,1)<1
    int_y(i,:) = int_y(i,:) + 1 - int_y(i,1);
end
if int_y(i,end)>d2
    int_y(i,:) = int_y(i,:) - (int_y(i,end)-d2);
end
Atemp = Atemp(int_x(i,:),int_y(i,:));
imagesc(int_x(i,:),int_y(i,:),Atemp); colormap(gca,'gray'); axis square; axis off; title(['ROI number ',num2str(i)])
end



