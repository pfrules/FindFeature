clc; clear; close all;

%files='C:\document\Findcircle\FindFeature\data\';
files='C:\document\Findcircle\FindFeature\datareal\';
rawImageFolder = [files '\figure\'];   % 输入文件路径

%  PixelPoints1=findpoint(rawImageFolder);%检测钢球坐标
%  save('point/PixelPoints1.mat','PixelPoints1')
load('point/PixelPoints1.mat')
PixelPoints2=CircleSubpixelDetection(rawImageFolder);
 save('point/PixelPoints2.mat','PixelPoints2')
load ([files '\PixelPoints.mat'])
PixelPoints3 = PhantomDetectIdealFunc();
upperCenterLast=[720,206];
lowerCenterLast=[720,1233];
num=size(PixelPoints1,1);
for i=1:num
[PixelPoints1(i).point,upperCenterLast,lowerCenterLast]=resort(PixelPoints1(i).point,upperCenterLast,lowerCenterLast);
end
upperCenterLast=[720,206];
lowerCenterLast=[720,1233];
for i=1:num
[PixelPoints2(i).point,upperCenterLast,lowerCenterLast]=resort(PixelPoints2(i).point,upperCenterLast,lowerCenterLast);
end
% figure;
% hold on
% for i=1:num
%     a(:,:)=PixelPoints(i).point;
%     b(:,:)=PixelPoints2(i).point;
%     plot(i,a(1,1),'r*')
%     plot(i,b(1,1),'b*')
% end
% for i=1:num
%     err(:,:,i)=PixelPoints(i).point-PixelPoints2(i).point;
% end

