clc; clear; close all;
format long;
filename='Data/';
pixel_scale = 1/(0.3);         % mm
upperCenterLast=[720,206];
lowerCenterLast=[720,1233];
flag=360;
% PixelPoints=get_pixel(pixel_scale,flag);
%yangyang
PixelPoints = PhantomDetectIdealFunc();
num=size(PixelPoints,1);
for i=1:num
%  [PixelPoints(i).point,upperCenterLast,lowerCenterLast]=resort(PixelPoints(i).point,upperCenterLast,lowerCenterLast);
PixelPoints(i).num=i;
end
% hold on;
% for i=1:720
% a=PixelPoints(i).point;
% plot(i,a(1,1),'r.');
% plot(i,a(1,2),'b.');
% end
FirstOutput=FirstCalc(PixelPoints,pixel_scale);%第一次参数计算
save([filename 'FirstOutput.mat'],'FirstOutput')
load([filename 'FirstOutput.mat']);


Parameter=SecondCalc(FirstOutput,pixel_scale);%转换到等中心点
save([filename 'Parameter.mat'],'Parameter')
load([filename 'Parameter.mat']);

%  hold on
% for i=1:size(Parameter,1)
% PsW=Parameter(i).PsW; 
% DW=Parameter(i).DW;
% plot3(PsW(1),PsW(2),PsW(3),'r.','MarkerSize',10)
% plot3(DW(1),DW(2),DW(3),'b.','MarkerSize',10)
% axis([xmin xmax ymin ymax -1 1]) 
% end
% for i=2:size(Parameter,1)
% plot(i,Parameter(i).t-Parameter(i-1).t,'r.')
% axis([0 720 0 1]) 
% end
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
t(i)=Parameter(i).t;
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
axis([0 360 0 1])
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
axis([0 360 0 1])  
hold off

% figure;
% hold on 
% plot(gantry,Pixel(gantry,1),'-r','LineWidth',2)
% plot(gantry,Pixel(gantry,2),'-b','LineWidth',2)
% legend('u','v');
% ylabel('偏移量（毫米）')
% xlabel('机架旋转角度（度）')
% axis([0 360 710 730])
% hold off
% %%
% figure;
% hold on
% plot(gantry,SID(gantry),'-r','LineWidth',2)
% plot(gantry,SAD(gantry),'-b','LineWidth',2)
% ylabel('源像距（毫米）')
% xlabel('机架旋转角度（度）')
% axis([0 360 2490 3510])  
% hold off
%%
% figure(6);
% hold on
% plot(gantry,SAD(gantry),'-r','LineWidth',2)
% ylabel('源轴距（毫米）')
% xlabel('机架旋转角度（度）')
% axis([0 360 2499.5 2500.5])  
% hold off
%%
% figure(7);
% hold on
% plot(gantry,t(gantry),'-r','LineWidth',2)
% ylabel('计算机架旋转角度（度）')
% xlabel('机架旋转角度（度）')
% axis([0 360 0 360])
% hold off

for i=1:num-1
t2(i)=Parameter(i+1).t-Parameter(i).t;
end
gantry=1:size(FirstOutput,1)-1;
figure;
hold on
plot(gantry,t2(gantry),'-r','LineWidth',2)
ylabel('扫描间隔角度（度）','FontSize',14)
xlabel('帧数','FontSize',14)
axis([0 360 0.95 1.05])  
hold off
