function hOutline = dispOutline(A_outline,d1,d2,disksize,col)
%adapted from epnev
%A_outline = full(A_or(:,k));
A_outline = reshape(A_outline,[d1,d2]);
A_outline(A_outline ~= 0) = 1;
A_outline = logical(A_outline);
se = strel('disk',disksize);
I = imopen(A_outline,se);
Ip = bwperim(I);
[yp,xp] = find(Ip == 1);

hOutline = plot(xp,yp,['.-',col],'MarkerSize',5);
end
