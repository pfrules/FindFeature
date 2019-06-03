function PixelPoints=newSubpixelDetection(rawImageFolder)
files = dir(strcat(rawImageFolder,'*.mat'));
LengthFiles = length(files);
PixelPoints = [];
% LengthFiles=1;
for p=1:LengthFiles
    file_image=strcat([rawImageFolder,files(p).name]);
    image=load(file_image);
    img=image.image;

    if size(files(1).name,2)>14
        [centers, radii,img_sub]  = ImagePreProcessor(img);
        img=65535-img;
        img = uint8(img / 256);
    else
        
        img= imrotate(img,90);
        [centers, radii] = imfindcircles(img,[8 12], 'ObjectPolarity', 'bright', ...
          'Sensitivity', 0.85, 'Method', 'twostage');
        img = uint8(img * 256);
    end
        time=3;
    [n,m]=size(img);
    [x,y]=meshgrid(1:m,1:n);%生成网格采样点
    [xi,yi]=meshgrid(1:time*m,1:time*n);

%     img=interp2(x,y,ff,xi/time,yi/time,'cubic');
    I=edge(interp2(x,y,double(img),xi/time,yi/time,'cubic'),'canny');
    figure;
    imshow(I);
    viscircles(time*centers,time*radii)
    figure;
    imshow(img);
    viscircles(centers,radii)
    BW=zeros(size(I));
    temp=ceil(1.414*max(time*radii));
    template=ones(2*temp,2*temp);
    for j=1:size(centers,1)
    BW(round(time*centers(j,2))-temp:round(time*centers(j,2))+temp-1,round(time*centers(j,1))-temp:round(time*centers(j,1))+temp-1)=template;
    end
    BW2=BW.*I;
% [lines,circles]=Edge_Segmentation(I);
% temp_img1=zeros(size(I));
% for k=1:length(lines)
%     edgePxiel=lines{k};
%     for i=1:size(edgePxiel,1)
%         temp_img1(edgePxiel(i,1),edgePxiel(i,2))=255;
%     end
% end
% figure,imshow(uint8(temp_img1))
% 
% temp_img2=zeros(size(I));
% for k=1:length(circles)
%     edgePxiel=circles{k};
%     for i=1:size(edgePxiel,1)
%         temp_img2(edgePxiel(i,1),edgePxiel(i,2))=255;
%     end
% end
% figure,imshow(uint8(temp_img2));
[centers, radii] = imfindcircles(I,[8 12], 'ObjectPolarity', 'bright', ...
          'Sensitivity', 0.85, 'Method', 'twostage');
viscircles(centers,radii)

 H=fspecial('gaussian',5,1);      % 高斯低通滤波器，用来过滤噪声
 img=imfilter(img,H);       %图像滤波、边界填充
%     figure;
%     imshow(img);
%     hold on
%     viscircles(centers,radii)
result=inmyway(img,centers, radii);
      

%     centers=round(centers);
     num=size(centers,1);
%     halfLength=ceil(sqrt(2)*max(radii));
%     for i=1:num
%     circle(:,:,i)=imcrop(img,[centers(i,:)-halfLength 2*halfLength 2*halfLength]);
%     end
%     centers=inMyWay(circle);
   centers=result; 
    
    n=size(img,2);
    if num<24
        continue;
    end
    xy1 =[];
    xy2 =[];
    for i=1:num
        c = centers(i,:);
        if(c(2)< n/2) %根据圆心的轴坐标分割成两组
            xy1 = [xy1;c];
        else
            xy2 = [xy2; c];
        end
    end

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
%     %显示编号和坐标
%     for i=1:size(xy1,1)
%         x1=round(xy1(i,1));
%         y1=round(xy1(i,2));
%         s1=sprintf(' %d -- (%d,%d)',i,x1,y1);
%         text(x1,y1,s1,'color','g','rotation',-45);
%         
%         x2=round(xy2(i,1));
%         y2=round(xy2(i,2));
%         s2=sprintf(' %d -- (%d,%d)',i,x2,y2);
%         text(x2,y2,s2,'color','g','rotation',-45)
%     end
    %将椭圆对应的坐标写入txt
    XY = [xy1, xy2];
  data.num=p-1;
  data.point=XY;
  PixelPoints=[PixelPoints;data];
        clc;
    disp(p);
end
end

function [seg_edges_line,seg_edges_circle]=Edge_Segmentation(input_img)
%函数功能：
%seg_edge=Edge_Segmentation(edge)
%实现边缘的分割
%函数输入：
%edge:边缘坐标
%函数输出：
%seg_edge:分割后的边缘段，为cell型
%seg_edge_num:分割后的边缘数:seg_edge_num(1)为水平边缘数，seg_edge_num(2)为垂直边缘数

if size(input_img,3)==3
    input_img=rbg2gray(input_img);
end

edge_img=input_img;
nrow=size(edge_img,1);
ncol=size(edge_img,2);

%进行霍夫变换
[H,T,R] = hough(edge_img,'Theta',-90:1:89);


%提取hough变换后的图像极值点并进行过滤
P=houghpeaks(H,10,'Threshold',0.2*max(H(:)));

%进行逐点扫描
scan_size=2;
seg_edges_line={};
for k=1:size(P,1)
    temp_edge=[];
    theta=T(P(k,2));
    rho=R(P(k,1));
    for x=1:ncol
        ideay=round((rho-x*cosd(theta))/sind(theta))+1
        if ideay<(size(edge_img,1)-scan_size) && ideay>scan_size
            patch=edge_img(ideay-scan_size:ideay+scan_size,x);
            realy=find(patch>0);
            if (length(realy)~=0)
                for ii=1:length(realy)
                    temp_edge(end+1,1:2)=[realy(ii)-scan_size-1+ideay,x];
                    edge_img(realy(ii)-scan_size-1+ideay,x)=0;
                end
            end
        end
    end
    if size(temp_edge,1)>15
        seg_edges_line{end+1}=temp_edge;
    end
    
end

temp_edges=seg_edges_line;
delete_len=25;
delete_gap=15;
for k=1:length(temp_edges);
    edge=temp_edges{k};
    real_start=1;
    real_len=[];
    len=1;
    for ii=2:size(edge,1)-1;
        if norm(edge(ii,:)-edge(ii+1,:))<delete_gap
            len=len+1;
        else
            if (len>delete_len)
                real_len(end+1)=len;
                real_start(end+1)=ii+1;
                len=1;
            else
                real_start(end)=ii+1;
                len=1;
            end
        end
    end
    if (len>delete_len)
        real_len(end+1)=len-1;    
    else
        real_start(end)=[];
    end
    
    temp_ed=[];
    for jj=1:length(real_start)
        temp_ed(end+1:end+real_len(jj),1:2)=edge(real_start(jj):real_start(jj)+real_len(jj)-1,:);
    end
    seg_edges_line{k}=temp_ed(1:end,:);
        
end

%======================================================================%
temp_img=input_img;
for k=1:length(seg_edges_line);
    edge=seg_edges_line{k};
    for i=1:size(edge,1)
        temp_img(edge(i,1),edge(i,2))=0;
    end
end

[label_img,num]=bwlabel(temp_img);
seg_edges_circle={};
for i=1:num
    edge_circle=[];
    [r,c]=find(label_img==i);
    if length(r)>delete_len
        edge_circle(end+1:end+length(r),1:2)=[r,c];
        seg_edges_circle{end+1}=edge_circle;
    end
end
end

function [centers, radii,img_sub] = ImagePreProcessor(fileLocation)
%%
%读取源图像
% File name example: "image_raw_00000.his"
I= fileLocation;
%figure;imshow(I);title('原图像');
%bw = roipoly(I);

%% 删除黑边
c1 = [90,60];
c2 = [1350,1400];
[m, n] = size(I);
pixelvalue = 65535;
imgcrop = I;
imgcrop(1:c1(1),:) = pixelvalue;
imgcrop(c2(1):m,:) = pixelvalue;
imgcrop(c1(1):c2(1),1:c1(2)) = pixelvalue;
imgcrop(c1(1):c2(1),c2(2):n) = pixelvalue;
imgcrop = pixelvalue -imgcrop;
%figure;imshow(imgcrop);title('删除黑边图像并反转图像');

%转换16bit图像为float
img = uint8(imgcrop / 256);
background = imopen(img,strel('disk',15));
img_sub = img-background;
%figure;imshow(img_sub);title('after');
bw = im2bw(img_sub,0.1); %设置合理的阈值保证获得质量好的二值化图像
bw1 = imopen(bw,strel('disk',5));
%figure;imshow(bw1);title('二值化图像');

%%
[centers, radii] = imfindcircles(bw1,[8 12], 'ObjectPolarity', 'bright', ...
          'Sensitivity', 0.85, 'Method', 'twostage'); %此处使用matlab自带函数进行圆检测，测试效果比较稳定

end