clear all;
close all;
x=-3;y=4;r=10;
t=0:0.1:360;
plot(t,-(x*cosd(t)+y*sind(t))+sqrt((x*cosd(t)+y*sind(t)).^2+(r^2-x^2-y^2)),'r');
hold on
phi=atand(x/y);
if y<0
    phi=phi+180;
end
if phi<0
    phi=360+phi;
end
plot(t,sqrt(5^2*sind(t+phi).^2+10^2-5^2)-5*sind(t+phi),'b');



%****************************************************************************************
%  
%                      ���������ź�Mix_Signal_1 ���ź� Mix_Signal_2 
%
%***************************************************************************************

Fs = 1000;                                                                        %������
N  = 1000;                                                                        %��������
n  = 0:N-1;
t   = 0:1/Fs:1-1/Fs;                                                            %ʱ������ 
Signal_Original_1 =sin(2*pi*10*t)+sin(2*pi*20*t)+sin(2*pi*30*t); 
Noise_White_1    = [0.3*randn(1,500), rand(1,500)];           %ǰ500���˹�ֲ�����������500����ȷֲ�������
Mix_Signal_1   = Signal_Original_1 + Noise_White_1;        %����Ļ���ź�

Signal_Original_2  =  [zeros(1,100), 20*ones(1,20), -2*ones(1,30), 5*ones(1,80), -5*ones(1,30), 9*ones(1,140), -4*ones(1,40), 3*ones(1,220), 12*ones(1,100), 5*ones(1,20), 25*ones(1,30), 7 *ones(1,190)]; 
Noise_White_2     =  0.5*randn(1,1000);                                 %��˹������
Mix_Signal_2        =  Signal_Original_2 + Noise_White_2;      %����Ļ���ź�

%****************************************************************************************
%  
%                �ź�Mix_Signal_1 �� Mix_Signal_2  �ֱ���������˹��ͨ�˲���butter�÷��μӸ�������
%
%***************************************************************************************

%����ź� Mix_Signal_1  ������˹��ͨ�˲�
figure(1);
Wc=2*50/Fs;                                          %��ֹƵ�� 50Hz
[b,a]=butter(4,Wc);
Signal_Filter=filter(b,a,Mix_Signal_1);

subplot(4,1,1);                                        %Mix_Signal_1 ԭʼ�ź�                 
plot(Mix_Signal_1);
axis([0,1000,-4,4]);
title('ԭʼ�ź� ');

subplot(4,1,2);                                        %Mix_Signal_1 ��ͨ�˲��˲����ź�  
plot(Signal_Filter);
axis([0,1000,-4,4]);
title('������˹��ͨ�˲����ź�');

%����ź� Mix_Signal_2  ������˹��ͨ�˲�
Wc=2*100/Fs;                                          %��ֹƵ�� 100Hz
[b,a]=butter(4,Wc);
Signal_Filter=filter(b,a,Mix_Signal_2);

subplot(4,1,3);                                        %Mix_Signal_2 ԭʼ�ź�                 
plot(Mix_Signal_2);
axis([0,1000,-10,30]);
title('ԭʼ�ź� ');

subplot(4,1,4);                                       %Mix_Signal_2 ��ͨ�˲��˲����ź�  
plot(Signal_Filter);
axis([0,1000,-10,30]);
title('������˹��ͨ�˲����ź�');

%****************************************************************************************
%  
%                �ź�Mix_Signal_1 �� Mix_Signal_2  �ֱ���FIR��ͨ�˲���
%
%***************************************************************************************

%����ź� Mix_Signal_1  FIR��ͨ�˲�
figure(2);
F   =  [0:0.05:0.95]; 
A  =  [1    1      0     0     0    0      0     0     0    0     0     0     0     0     0     0    0   0   0   0] ;
b  =  firls(20,F,A);
Signal_Filter = filter(b,1,Mix_Signal_1);

subplot(4,1,1);                                          %Mix_Signal_1 ԭʼ�ź�                 
plot(Mix_Signal_1);
axis([0,1000,-4,4]);
title('ԭʼ�ź� ');

subplot(4,1,2);                                          %Mix_Signal_1 FIR��ͨ�˲��˲����ź�  
plot(Signal_Filter);
axis([0,1000,-5,5]);
title('FIR��ͨ�˲�����ź�');

%����ź� Mix_Signal_2  FIR��ͨ�˲�
F   =  [0:0.05:0.95]; 
A  =  [1    1      1     1     1    0      0    0     0    0     0     0     0     0     0     0    0   0   0   0] ;
b  =  firls(20,F,A);
Signal_Filter = filter(b,1,Mix_Signal_2);
subplot(4,1,3);                                          %Mix_Signal_2 ԭʼ�ź�                 
plot(Mix_Signal_2);
axis([0,1000,-10,30]);
title('ԭʼ�ź� ');

subplot(4,1,4);                                          %Mix_Signal_2 FIR��ͨ�˲��˲����ź�  
plot(Signal_Filter);
axis([0,1000,-10,30]);
title('FIR��ͨ�˲�����ź�');

%****************************************************************************************
%  
%                �ź�Mix_Signal_1 �� Mix_Signal_2  �ֱ����ƶ�ƽ���˲�
%
%***************************************************************************************

%����ź� Mix_Signal_1  �ƶ�ƽ���˲�
figure(3);
b  =  [1 1 1 1 1 1]/6;
Signal_Filter = filter(b,1,Mix_Signal_1);

subplot(4,1,1);                                          %Mix_Signal_1 ԭʼ�ź�                 
plot(Mix_Signal_1);
axis([0,1000,-4,4]);
title('ԭʼ�ź� ');

subplot(4,1,2);                                          %Mix_Signal_1 �ƶ�ƽ���˲����ź�  
plot(Signal_Filter);
axis([0,1000,-4,4]);
title('�ƶ�ƽ���˲�����ź�');

%����ź� Mix_Signal_2  �ƶ�ƽ���˲�
b  =  [1 1 1 1 1 1]/6;
Signal_Filter = filter(b,1,Mix_Signal_2);
subplot(4,1,3);                                          %Mix_Signal_2 ԭʼ�ź�                 
plot(Mix_Signal_2);
axis([0,1000,-10,30]);
title('ԭʼ�ź� ');

subplot(4,1,4);                                          %Mix_Signal_2 �ƶ�ƽ���˲����ź�  
plot(Signal_Filter);
axis([0,1000,-10,30]);
title('�ƶ�ƽ���˲�����ź�');

%****************************************************************************************
%  
%                �ź�Mix_Signal_1 �� Mix_Signal_2  �ֱ�����ֵ�˲�
%
%***************************************************************************************

%����ź� Mix_Signal_1  ��ֵ�˲�
figure(4);
Signal_Filter=medfilt1(Mix_Signal_1,10);

subplot(4,1,1);                                          %Mix_Signal_1 ԭʼ�ź�                 
plot(Mix_Signal_1);
axis([0,1000,-5,5]);
title('ԭʼ�ź� ');

subplot(4,1,2);                                          %Mix_Signal_1 ��ֵ�˲����ź�  
plot(Signal_Filter);
axis([0,1000,-5,5]);
title('��ֵ�˲�����ź�');

%����ź� Mix_Signal_2  ��ֵ�˲�
Signal_Filter=medfilt1(Mix_Signal_2,10);
subplot(4,1,3);                                          %Mix_Signal_2 ԭʼ�ź�                 
plot(Mix_Signal_2);
axis([0,1000,-10,30]);
title('ԭʼ�ź� ');

subplot(4,1,4);                                          %Mix_Signal_2 ��ֵ�˲����ź�  
plot(Signal_Filter);
axis([0,1000,-10,30]);
title('��ֵ�˲�����ź�');

%****************************************************************************************
%  
%                �ź�Mix_Signal_1 �� Mix_Signal_2  �ֱ���ά���˲�
%
%***************************************************************************************

%����ź� Mix_Signal_1  ά���˲�
figure(5);
Rxx=xcorr(Mix_Signal_1,Mix_Signal_1);              %�õ�����źŵ�����غ���
M=100;                                                             %ά���˲�������
for i=1:M                                                           %�õ�����źŵ�����ؾ���
    for j=1:M
        rxx(i,j)=Rxx(abs(j-i)+N);
    end
end
Rxy=xcorr(Mix_Signal_1,Signal_Original_1);       %�õ�����źź�ԭ�źŵĻ���غ���
for i=1:M
    rxy(i)=Rxy(i+N-1);
end                                                                  %�õ�����źź�ԭ�źŵĻ��������
h = inv(rxx)*rxy';                                               %�õ���Ҫ�漰��wiener�˲���ϵ��
Signal_Filter=filter(h,1, Mix_Signal_1);               %�������ź�ͨ��ά���˲���

subplot(4,1,1);                                                   %Mix_Signal_1 ԭʼ�ź�                 
plot(Mix_Signal_1);
axis([0,1000,-5,5]);
title('ԭʼ�ź� ');

subplot(4,1,2);                                                   %Mix_Signal_1 ά���˲����ź�  
plot(Signal_Filter);
axis([0,1000,-5,5]);
title('ά���˲�����ź�');

%����ź� Mix_Signal_2  ά���˲�
Rxx=xcorr(Mix_Signal_2,Mix_Signal_2);              %�õ�����źŵ�����غ���
M=500;                                                             %ά���˲�������
for i=1:M                                                           %�õ�����źŵ�����ؾ���
    for j=1:M
        rxx(i,j)=Rxx(abs(j-i)+N);
    end
end
Rxy=xcorr(Mix_Signal_2,Signal_Original_2);       %�õ�����źź�ԭ�źŵĻ���غ���
for i=1:M
    rxy(i)=Rxy(i+N-1);
end                                                                  %�õ�����źź�ԭ�źŵĻ��������
h=inv(rxx)*rxy';                                               %�õ���Ҫ�漰��wiener�˲���ϵ��
Signal_Filter=filter(h,1, Mix_Signal_2);             %�������ź�ͨ��ά���˲���

subplot(4,1,3);                                                  %Mix_Signal_2 ԭʼ�ź�                 
plot(Mix_Signal_2);
axis([0,1000,-10,30]);
title('ԭʼ�ź� ');

subplot(4,1,4);                                                  %Mix_Signal_2 ά���˲����ź�  
plot(Signal_Filter);
axis([0,1000,-10,30]);
title('ά���˲�����ź�');

%****************************************************************************************
%  
%                �ź�Mix_Signal_1 �� Mix_Signal_2  �ֱ�������Ӧ�˲�
%
%***************************************************************************************

%����ź� Mix_Signal_1 ����Ӧ�˲�
figure(6);
N=1000;                                             %�����źų�������N
k=100;                                                  %ʱ���ͷLMS�㷨�˲�������
u=0.001;                                             %��������

%���ó�ֵ
yn_1=zeros(1,N);                                  %output signal
yn_1(1:k)=Mix_Signal_1(1:k);                 %�������ź�SignalAddNoise��ǰk��ֵ��Ϊ���yn_1��ǰk��ֵ
w=zeros(1,k);                                        %���ó�ͷ��Ȩ��ֵ
e=zeros(1,N);                                        %����ź�

%��LMS�㷨�����˲�
for i=(k+1):N
        XN=Mix_Signal_1((i-k+1):(i));
        yn_1(i)=w*XN';
        e(i)=Signal_Original_1(i)-yn_1(i);
        w=w+2*u*e(i)*XN;
end

subplot(4,1,1);
plot(Mix_Signal_1);                               %Mix_Signal_1 ԭʼ�ź�
axis([k+1,1000,-4,4]);
title('ԭʼ�ź�');

subplot(4,1,2);
plot(yn_1);                                            %Mix_Signal_1 ����Ӧ�˲����ź�
axis([k+1,1000,-4,4]);
title('����Ӧ�˲����ź�');

%����ź� Mix_Signal_2 ����Ӧ�˲�
N=1000;                                             %�����źų�������N
k=500;                                                %ʱ���ͷLMS�㷨�˲�������
u=0.000011;                                        %��������

%���ó�ֵ
yn_1=zeros(1,N);                                   %output signal
yn_1(1:k)=Mix_Signal_2(1:k);                  %�������ź�SignalAddNoise��ǰk��ֵ��Ϊ���yn_1��ǰk��ֵ
w=zeros(1,k);                                        %���ó�ͷ��Ȩ��ֵ
e=zeros(1,N);                                        %����ź�

%��LMS�㷨�����˲�
for i=(k+1):N
        XN=Mix_Signal_2((i-k+1):(i));
        yn_1(i)=w*XN';
        e(i)=Signal_Original_2(i)-yn_1(i);
        w=w+2*u*e(i)*XN;
end

subplot(4,1,3);
plot(Mix_Signal_2);                               %Mix_Signal_1 ԭʼ�ź�
axis([k+1,1000,-10,30]);
title('ԭʼ�ź�');

subplot(4,1,4);
plot(yn_1);                                            %Mix_Signal_1 ����Ӧ�˲����ź�
axis([k+1,1000,-10,30]);
title('����Ӧ�˲����ź�');

%****************************************************************************************
%  
%                �ź�Mix_Signal_1 �� Mix_Signal_2  �ֱ���С���˲�
%
%***************************************************************************************

%����ź� Mix_Signal_1  С���˲�
figure(7);
subplot(4,1,1);
plot(Mix_Signal_1);                                 %Mix_Signal_1 ԭʼ�ź�
axis([0,1000,-5,5]);
title('ԭʼ�ź� ');

subplot(4,1,2);
[xd,cxd,lxd] = wden(Mix_Signal_1,'sqtwolog','s','one',2,'db3');
plot(xd);                                                 %Mix_Signal_1 С���˲����ź�
axis([0,1000,-5,5]);
title('С���˲����ź� ');

%����ź� Mix_Signal_2  С���˲�
subplot(4,1,3);
plot(Mix_Signal_2);                                 %Mix_Signal_2 ԭʼ�ź�
axis([0,1000,-10,30]);
title('ԭʼ�ź� ');

subplot(4,1,4);
[xd,cxd,lxd] = wden(Mix_Signal_2,'sqtwolog','h','sln',3,'db3');
plot(xd);                                                %Mix_Signal_2 С���˲����ź�
axis([0,1000,-10,30]);
title('С���˲����ź� ');