function PixelPoints=subpixelDetection(rawImageFolder)
files = dir(strcat(rawImageFolder,'*.mat'));
LengthFiles = length(files);
PixelPoints = [];
% LengthFiles=1;
for p=1:LengthFiles
    file_image=strcat([rawImageFolder,files(p).name]);
    image=load(file_image);
    if size(files(1).name,2)>14
        [centers, radii,img_sub]  = ImagePreProcessor(image);
        img=65535-image.image;
        img = uint8(img / 256);
    else
        
        img= imrotate(image.image,90);
        [centers, radii] = imfindcircles(img,[8 12], 'ObjectPolarity', 'bright', ...
          'Sensitivity', 0.85, 'Method', 'twostage');
        img = uint8(img * 256);
    end

%     figure;
%     imshow(img);
%     hold on
inmyway(img,centers, radii);
      
%     viscircles(centers,radii)
    centers=round(centers);
    num=size(centers,1);
    halfLength=ceil(sqrt(2)*max(radii));
    for i=1:num
    circle(:,:,i)=imcrop(img,[centers(i,:)-halfLength 2*halfLength 2*halfLength]);
    end
    centers=inMyWay(circle);
    
    
    n=size(img,2);
    if num<24
        continue;
    end
    xy1 =[];
    xy2 =[];
    for i=1:num
        c = centers(i,:);
        if(c(2)< n/2) %����Բ�ĵ�������ָ������
            xy1 = [xy1;c];
        else
            xy2 = [xy2; c];
        end
    end

    center_x1 = xy1(:,1);
    center_y1 = xy1(:,2);
    center_x2 = xy2(:,1);
    center_y2 = xy2(:,2);
    k1 = convhull(center_x1, center_y1); %ע��˴���13��ֵ
    k2 = convhull(center_x2, center_y2);
    k2 = k2(end:-1:1);
    xy1 = xy1(k1(1:end-1),:);
    xy2 = xy2(k2(1:end-1),:);
    [min1,index1] = min(xy1);
    [min2,index2] = min(xy2);

    %��������
    xy1 = [xy1(index1:end, :); xy1(1:index1-1, :)];
    xy2 = [xy2(index2:end, :); xy2(1:index2-1, :)];
%     %��ʾ��ź�����
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
    %����Բ��Ӧ������д��txt
    XY = [xy1, xy2];
  data.num=p-1;
  data.point=XY;
  PixelPoints=[PixelPoints;data];
        clc;
    disp(p);
end
end

function centers=inMyWay(circle)
for i=1:size(circle,3)
BW=edge(circle(:,:,i),'canny');

%figure(i);m=1:35;n=1:35;meshc(m,n,circle(m,n,i))
% [oimg,X,Y]=BaseZernike1(circle(:,:,i));
%O=Edgedectect(circle(:,:,i));
[oimg,X,Y]=myZernike7(circle(:,:,i));
end
end
function [centers, radii,img_sub] = ImagePreProcessor(fileLocation)
%%
%��ȡԴͼ��
% File name example: "image_raw_00000.his"
I= fileLocation.image;
%figure;imshow(I);title('ԭͼ��');
%bw = roipoly(I);

%% ɾ���ڱ�
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
%figure;imshow(imgcrop);title('ɾ���ڱ�ͼ�񲢷�תͼ��');

%ת��16bitͼ��Ϊfloat
img = uint8(imgcrop / 256);
background = imopen(img,strel('disk',15));
img_sub = img-background;
%figure;imshow(img_sub);title('after');
bw = im2bw(img_sub,0.1); %���ú������ֵ��֤��������õĶ�ֵ��ͼ��
bw1 = imopen(bw,strel('disk',5));
%figure;imshow(bw1);title('��ֵ��ͼ��');

%%
[centers, radii] = imfindcircles(bw1,[8 12], 'ObjectPolarity', 'bright', ...
          'Sensitivity', 0.85, 'Method', 'twostage'); %�˴�ʹ��matlab�Դ���������Բ��⣬����Ч���Ƚ��ȶ�

end