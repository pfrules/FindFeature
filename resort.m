function [XY,upperCenterCurrent,lowerCenterCurrent]=resort(centers,upperCenterLast,lowerCenterLast)


%     num=size(centers,1);
%     xy1 =centers(centers(:,2)<mean(centers(:,2)),:);
%     xy2 =centers(centers(:,2)>mean(centers(:,2)),:);
    xy1 =centers(:,1:2);
    xy2 =centers(:,3:4);

    center_x1 = xy1(:,1);
    center_y1 = xy1(:,2);
    center_x2 = xy2(:,1);
    center_y2 = xy2(:,2);
    k1 = convhull(center_x1, center_y1); %注意此处有13个值
    k2 = convhull(center_x2, center_y2);
    k2 = k2(end:-1:1);
    xy1 = xy1(k1(1:end-1),:);
    xy2 = xy2(k2(1:end-1),:);

[min1,index1] = min(xy1);
[min2,index2] = min(xy2);

%重新排列
xy1 = [xy1(index1:end, :); xy1(1:index1-1, :)];
xy2 = [xy2(index2:end, :); xy2(1:index2-1, :)];
[upperCenterCurrent,index1] = NearestObject(xy1,upperCenterLast);
[lowerCenterCurrent,index2] = NearestObject(xy2,lowerCenterLast);

%重新排列
xy1 = [xy1(index1:end, :); xy1(1:index1-1, :)];
xy2 = [xy2(index2:end, :); xy2(1:index2-1, :)];

%将椭圆对应的坐标写入txt
XY = [xy1, xy2];
end
    
function [posi, index] = NearestObject(XY, p)
%
%Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%        p 为指定的搜索起始坐标
%Output: p = [x, y] 为最近的钢球的位置
%        index = i 为对应的XY的下标

[m,n] = size(XY);
dis=sum((repmat(p,m,1)-XY).^2,2);
[index,~]=find(dis==min(dis));
posi = XY(index, :);
end