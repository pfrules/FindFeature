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
    boxNum=size(centers,1);
    %插值倍数
        time=2;
    [n,m]=size(img);
    [x,y]=meshgrid(1:m,1:n);%生成网格采样点
    [xi,yi]=meshgrid(1:time*m,1:time*n);

    interpImg=uint8(interp2(x,y,double(img),xi/time,yi/time,'cubic'));
    I=edge(interpImg,'canny');%,linear
%     figure;
%     imshow(I);
%     viscircles(time*centers,time*radii)
%     figure;
%     imshow(img);
%     viscircles(centers,radii)
    
    
    BW=zeros(size(I));
    temp=ceil(1.414*max(time*radii));
    template=ones(2*temp,2*temp);
    for j=1:boxNum
    BW(round(time*centers(j,2))-temp:round(time*centers(j,2))+temp-1,round(time*centers(j,1))-temp:round(time*centers(j,1))+temp-1)=template;
    end
    BW2=uint8(BW.*I);
    
%     显示第一次处理的边缘信息
%     figure;
%     imshow(BW2);
%     viscircles(time*centers,time*radii)
    
    edgeCircle=EdgeSegmentation(BW2,boxNum);
BW=zeros(size(I));
figure;
imshow(interpImg);
viscircles(time*centers,time*radii)
hold on
for i=1:boxNum
    fitcircle=[];
    boxEdge=edgeCircle(i).PixelList;
    boxEdgeNum=edgeCircle(i).Area;
    for j=1:boxEdgeNum
        BW(boxEdge(j,2),boxEdge(j,1))=1;
        %plot(boxEdge(j,1),boxEdge(j,2),'y*')
    end
    %%直接圆拟合
    [par,boxEdge]=excludepoint(boxEdge,time);
%     h2=ezplot(@(x,y)x^2+y^2+par(2)*x+par(3)*y+par(4),[min(boxEdge(:,1))-3,max(boxEdge(:,1))+3,min(boxEdge(:,2))-3,max(boxEdge(:,2))+3])
%     set(h2,'Color','k');

    %%在圆内取一点为中心，在360度方向上计算距离的曲线
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
    [Re, center, vertex]= MyEllipseDirectFit(boxEdge);
    eq0=@(x,y) Re(1)*x^2 + Re(2)*x*y + Re(3)*y^2 +Re(4)*x + Re(5)*y + Re(6);
    h=ezplot(eq0,[min(boxEdge(:,1))-3,max(boxEdge(:,1))+3,min(boxEdge(:,2))-3,max(boxEdge(:,2))+3]);
    set(h,'Color','b');
end


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

function [par,boxEdge]=excludepoint(boxEdge,time)
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

function edgeCircle=EdgeSegmentation(inputImg,numbox)
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

edgeCircle=img_reg(1:numbox);
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