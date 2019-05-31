clc; clear; close all;

%files='C:\document\Findcircle\FindFeature\data\';
files='C:\document\Findcircle\FindFeature\datareal\';
rawImageFolder = [files '\figure\'];   % �����ļ�·��

%  PixelPoints1=findpoint(rawImageFolder);%����������
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
FirstOutput=FirstCalc(PixelPoints,pixel_scale);%��һ�β�������
save([files 'FirstOutput.mat'],'FirstOutput')
load([files 'FirstOutput.mat']);


Parameter=SecondCalc(FirstOutput,pixel_scale);%ת���������ĵ�
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
ylabel('ƽ����ת�Ƕȣ��ȣ�','FontSize',14)
xlabel('֡��','FontSize',14)
axis([0 360 -1 1])
hold off

% figure(2);
% hold on 
% plot(gantry,set(gantry,1),'-r','LineWidth',2)
% plot(gantry,set(gantry,2),'-b','LineWidth',2)
% legend('u','v');
% ylabel('ƫ���������ף�')
% xlabel('������ת�Ƕȣ��ȣ�')
% axis([0 360 -10 10])
% hold off



%%
figure;
hold on
plot(gantry,theta2(gantry),'-r','LineWidth',2)
plot(gantry,phi2(gantry),'-b','LineWidth',2)
plot(gantry,eta2(gantry),'-k','LineWidth',2)
legend('theta','phi','eta');
ylabel('ƽ����ת�Ƕȣ��ȣ�','FontSize',14)
xlabel('֡��','FontSize',14)
axis([0 360 -1 1])  
hold off
