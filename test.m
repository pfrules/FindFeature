clc; clear; close all;

%files='C:\document\Findcircle\FindFeature\data\';
files='C:\document\Findcircle\FindFeature\datareal\';
rawImageFolder = [files '\figure\'];   % 输入文件路径

%  PixelPoints1=findpoint(rawImageFolder);%检测钢球坐标
%  save('point/PixelPoints1.mat','PixelPoints1')
load('point/PixelPoints1.mat')
 PixelPoints2=subpixelDetection(rawImageFolder);
 save('point/PixelPoints2.mat','PixelPoints2')
load ([files '\PixelPoints.mat'])
PixelPoints3 = PhantomDetectIdealFunc();
upperCenterLast=[720,206];
lowerCenterLast=[720,1233];
num=size(PixelPoints,1);
for i=1:num
[PixelPoints(i).point,upperCenterLast,lowerCenterLast]=resort(PixelPoints(i).point,upperCenterLast,lowerCenterLast);
PixelPoints(i).num=i;
end
% upperCenterLast=[720,206];
% lowerCenterLast=[720,1233];
% for i=1:num
% [PixelPoints2(i).point,upperCenterLast,lowerCenterLast]=resort(PixelPoints2(i).point,upperCenterLast,lowerCenterLast);
% PixelPoints2(i).num=i;
% end

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
pixel_scale = 1/(0.3); 
FirstOutput=FirstCalc(PixelPoints,pixel_scale);%第一次参数计算
save([files 'FirstOutput.mat'],'FirstOutput')
load([files 'FirstOutput.mat']);


Parameter=SecondCalc(FirstOutput,pixel_scale);%转换到等中心点
save([files 'Parameter.mat'],'Parameter')
load([files 'Parameter.mat']);
num=size(FirstOutput,1);
for i=1:num
theta1(i)=FirstOutput(i).theta;
phi1(i)=FirstOutput(i).phi;
eta1(i)=FirstOutput(i).eta;

end
for i=1:num
theta2(i)=Parameter(i).theta;
phi2(i)=Parameter(i).phi;
eta2(i)=Parameter(i).eta;
SID(i)=Parameter(i).SID;
end
gantry=1:num;
figure(1);
hold on 
plot(gantry,theta1(gantry),'-r','LineWidth',2)
plot(gantry,phi1(gantry),'-b','LineWidth',2)
plot(gantry,eta1(gantry),'-k','LineWidth',2)
legend('theta','phi','eta');
ylabel('平板旋转角度（度）','FontSize',14)
xlabel('帧数','FontSize',14)
axis([0 360 -1 1])
hold off

% figure(2);
% hold on 
% plot(gantry,set(gantry,1),'-r','LineWidth',2)
% plot(gantry,set(gantry,2),'-b','LineWidth',2)
% legend('u','v');
% ylabel('偏移量（毫米）')
% xlabel('机架旋转角度（度）')
% axis([0 360 -10 10])
% hold off



%%
figure;
hold on
plot(gantry,theta2(gantry),'-r','LineWidth',2)
plot(gantry,phi2(gantry),'-b','LineWidth',2)
plot(gantry,eta2(gantry),'-k','LineWidth',2)
legend('theta','phi','eta');
ylabel('平板旋转角度（度）','FontSize',14)
xlabel('帧数','FontSize',14)
axis([0 360 -1 1])  
hold off
