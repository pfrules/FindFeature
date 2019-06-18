function PixelPoints=CircleSubpixelDetection(rawImageFolder)
files = dir(strcat(rawImageFolder,'*.mat'));
LengthFiles = length(files);
PixelPoints = [];

show=0;
for p=1:LengthFiles
    file_image=strcat([rawImageFolder,files(p).name]);
    image=load(file_image);
    img=image.image;

    if size(files(1).name,2)>14
        [centers, radii]  = ImagePreProcessor(img);
        img=65535-img;
        img = uint8(img / 256);
    else
        img= imrotate(img,90);
        [centers, radii] = imfindcircles(img,[8 12], 'ObjectPolarity', 'bright', ...
          'Sensitivity', 0.85, 'Method', 'twostage');
        img = uint8(img * 256);
    end
   
% result=lowway(img,centers, radii); %lowway使用像素拟合的方法，实验证明效果不好
    boxNum=size(centers,1);
%%
%图像插值
     multiple=1;  %图像插值倍数
     
    [n,m]=size(img);
    [x,y]=meshgrid(1:m,1:n);%生成网格采样点
    [xi,yi]=meshgrid(1:multiple*m,1:multiple*n);
    interpImg=uint8(interp2(x,y,double(img),xi/multiple,yi/multiple,'cubic'));%cubic双立方插值后的结果%,linear
%%%%%%%%%%%%%%%%
%%zernike矩检测亚像素，返回亚像素阶跃坐标
%     H=fspecial('gaussian',5,1);      % 高斯低通滤波器，用来过滤噪声
%     interpImg2=imfilter(interpImg,H);       %图像滤波、边界填充
%     [col,row,Zernikek,Zernikel]=BaseZernike(interpImg2);
%     [sizerow,sizecol]=size(col);
%%%%%%%%%%%%%%%

    I=edge(interpImg,'canny');
%     figure;
%     imshow(I);
%     viscircles(time*centers,time*radii)
%     figure;
%     imshow(img);
%     viscircles(centers,radii)
    
%%
%删除其他部分干扰项
    BW=zeros(size(I));
    temp=ceil(1.414*max(multiple*radii));
    template=ones(2*temp,2*temp);
    for j=1:boxNum
    BW(round(multiple*centers(j,2))-temp:round(multiple*centers(j,2))+temp-1,round(multiple*centers(j,1))-temp:round(multiple*centers(j,1))+temp-1)=template;
    end
    BW2=uint8(BW.*I);
    
%     显示第一次处理的边缘信息
if show
 figure;
imshow(interpImg);
viscircles(multiple*centers,multiple*radii)
hold on  
end
%%
%粗略分割钢球边缘像素值
edgeCircle=EdgeSegmentation(BW2,centers, radii,boxNum,multiple);
BW=zeros(size(I));

for i=1:boxNum
    boxEdge=edgeCircle(i).PixelList;
    boxEdgeNum=edgeCircle(i).Area;
    for j=1:boxEdgeNum
        BW(boxEdge(j,2),boxEdge(j,1))=1;
        if show
         plot(boxEdge(j,1),boxEdge(j,2),'y*')           
        end

    end
%%
%删除每个钢球的干扰像素点的两种方法
    %%1、直接圆拟合
    [par,boxEdge]=excludepoint(boxEdge,multiple);
%     h2=ezplot(@(x,y)x^2+y^2+par(2)*x+par(3)*y+par(4),[min(boxEdge(:,1))-3,max(boxEdge(:,1))+3,min(boxEdge(:,2))-3,max(boxEdge(:,2))+3])
%     set(h2,'Color','k');

    %%2、在圆内取一点为中心，在360度方向上计算距离的曲线
%     blackBoxEdge=boxEdge-repmat(mean(boxEdge),size(boxEdge,1),1);
%     angle=atand(blackBoxEdge(:,2)./blackBoxEdge(:,1));
%     angle(blackBoxEdge(:,1)<0)=180+angle(blackBoxEdge(:,1)<0);
%     angle(angle<0)=360+angle(angle<0);
%     rol=sqrt(blackBoxEdge(:,1).^2+blackBoxEdge(:,2).^2);
%     Signal=[angle,rol,[1:size(boxEdge)]'];
%     Signal=sortrows(Signal,1);
% %     figure(i+1)
% %     plot(Signal(:,1),Signal(:,2),'r*');
%     [parameter,Signal]=optimization(Signal,mean(time*radii),time);
%     hold on
%     x0=parameter(1);
%     y0=parameter(2);
%     r=parameter(3);
%     phi=atand(x0/y0);
% if y0<0
%     phi=phi+180;
% end
% if phi<0
%     phi=360+phi;
% end
% %     plot(Signal(:,1),sqrt((x0^2+y0^2)*sind(Signal(:,1)+phi).^2+r^2-(x0^2+y0^2))-sqrt(x0^2+y0^2)*sind(Signal(:,1)+phi),'-b');
% %     plot(Signal(:,1),Signal(:,2),'o');
%     boxEdge=boxEdge(Signal(:,3),:);

    %plot(boxEdge(:,1),boxEdge(:,2),'ok');
  
%zernike矩，取亚像素坐标，效果并不好
%     cols=col((boxEdge(:,2)-1)*sizecol+boxEdge(:,1));
%     rows=row((boxEdge(:,2)-1)*sizecol+boxEdge(:,1));
%     boxE2dgek=Zernikek((boxEdge(:,2)-1)*sizecol+boxEdge(:,1));
%     boxEdgel=Zernikel((boxEdge(:,2)-1)*sizecol+boxEdge(:,1));
%     boxEdge=[rows,cols];
%     boxEdge3=round(boxEdge2);
%     plot(boxEdge(:,1),boxEdge(:,2),'r*');
%     plot(boxEdge3(:,1),boxEdge3(:,2),'oy');
    [Re, center, vertex]= MyEllipseDirectFit(boxEdge);%椭圆拟合
    if show
    eq0=@(x,y) Re(1)*x^2 + Re(2)*x*y + Re(3)*y^2 +Re(4)*x + Re(5)*y + Re(6);
    h=ezplot(eq0,[min(boxEdge(:,1))-3,max(boxEdge(:,1))+3,min(boxEdge(:,2))-3,max(boxEdge(:,2))+3]);
    set(h,'Color','b');     
    end
result(i,:)=center;
end

    if show
        plot(result(:,1),result(:,2),'r*');
    end
   centers=result/multiple; 
      num=size(centers,1);   
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

function [par,boxEdge]=excludepoint(boxEdge,time)
%圆拟合，删除圆外像素点多余的
    fitcircle(:,1)=boxEdge(:,1).^2+boxEdge(:,2).^2;
    fitcircle(:,2)=boxEdge(:,1);
    fitcircle(:,3)=boxEdge(:,2);
    fitcircle(:,4)=1;
    [U,D,V]=svd(fitcircle);
    par=V(:,end)/V(1,end);
    
    x0=-par(2)/2;
    y0=-par(3)/2;
    r=sqrt(par(2)^2+par(3)^2-4*par(4))/2;
    dis=sqrt(sum((boxEdge-repmat([x0,y0],size(boxEdge,1),1)).^2,2))-r;
    if ~isempty(find(abs(dis)>log2(time)+1.5))
        boxEdge(abs(dis)==max(abs(dis)),:)=[];
        [par,boxEdge]=excludepoint(boxEdge,time);
    end
end


function [parameter,Signal]=optimization(Signal,r,time)
%LM非线性拟合
x0=0.1;y0=0.1;
    para=[x0,y0,r];
 %    options = optimset('LargeScale','off','LevenbergMarquardt','on');
      options = optimset('Display','off','TolX',eps,'TolFun',eps,'LargeScale','off','Algorithm','Levenberg-Marquardt');
    [parameter,resnorm,residual,exitflag,output]  = lsqnonlin( @myop, para, [],[],options, Signal);
% display the result
x0=parameter(1);
y0=parameter(2);
r=parameter(3);
phi=atand(x0/y0);
if y0<0
    phi=phi+180;
end
if phi<0
    phi=360+phi;
end
%%
f=sqrt((x0^2+y0^2)*sind(Signal(:,1)+phi).^2+r^2-(x0^2+y0^2))-sqrt(x0^2+y0^2)*sind(Signal(:,1)+phi)-Signal(:,2);
if ~isempty(find(abs(f)>log2(time)+1.5))
    Signal(abs(f)==max(abs(f)),:)=[];
    [parameter,Signal]=optimization(Signal,r,time);
end
end


function f = myop(params, Signal)
x0=params(1);
y0=params(2);
r=params(3);
phi=atand(x0/y0);
if y0<0
    phi=phi+180;
end
if phi<0
    phi=360+phi;
end
%%
f=sqrt((x0^2+y0^2)*sind(Signal(:,1)+phi).^2+r^2-(x0^2+y0^2))-sqrt(x0^2+y0^2)*sind(Signal(:,1)+phi)-Signal(:,2);
end

function edgeCircle=EdgeSegmentation(inputImg,centers, radii,numbox,multiple)
%函数功能：
%提取圆的部分

if size(inputImg,3)==3
    inputImg=rbg2gray(inputImg);
end
[label_img,num]=bwlabel(inputImg);
img_reg = regionprops(label_img, 'Area','PixelList'); 
% for i=1:num
%     [x,y]=find(label_img==i);
%     area(i,1)=i;
%     area(i,2)=length(x);
% end
% for i=1:num
%     img_reg(i).num=i;
% end
for i=1:num
    for j=i+1:num
        if img_reg(i).Area<img_reg(j).Area
            temp=img_reg(i);
            img_reg(i)=img_reg(j);
            img_reg(j)=temp;
        end
    end
end
edge0=[];
for i=1:numbox
    edge0=[edge0;img_reg(i).PixelList];
end
for i=1:numbox
    dis=abs(sqrt(sum((edge0-repmat(multiple*centers(i,:),size(edge0,1),1)).^2,2))-multiple*radii(i));
    edgeCircle(i).PixelList=edge0(dis<multiple*2,:);
    edgeCircle(i).Area=size(edgeCircle(i).PixelList,1);
end
% plot(edge0(:,1),edge0(:,2),'y*')
end

function [centers, radii] = ImagePreProcessor(fileLocation)
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
[label_img,num]=bwlabel(bw1);
img_reg = regionprops(label_img, 'Area','Centroid'); 
for i=1:num
    for j=i+1:num
        if img_reg(i).Area<img_reg(j).Area
            temp=img_reg(i);
            img_reg(i)=img_reg(j);
            img_reg(j)=temp;
        end
    end
end
if num>24
    img_reg=img_reg(1:24);
end
% for i=1:size(img_reg,1)
%     radii(i,:)=sqrt(img_reg(i).Area/pi);
%     centers(i,:)=img_reg(i).Centroid;
% end
[centers, radii] = imfindcircles(bw1,[8 12], 'ObjectPolarity', 'bright', ...
          'Sensitivity', 0.85, 'Method', 'twostage'); %此处使用matlab自带函数进行圆检测，测试效果比较稳定
point=[centers,radii];
%部分图像会出现一个钢球检测出两个圆的情况，合并多余的
if size(point,1)>24
    point=deletedata(point);
end
centers=point(:,1:2);
radii=point(:,3);
end

function point=deletedata(point)
num=size(point,1);
for i=1:num
    dis=sqrt(sum((point(:,1:2)-repmat(point(i,1:2),num,1)).^2,2));
    sorted=sort(dis);
    point(i,4)=sorted(2);
end
point=sortrows(point,4);
point(2,:)=mean(point(1:2,:));
point=point(2:end,1:3);
if size(point,1)>24
    point=deletedata(point);
end
end