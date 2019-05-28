function inmyway(img,centers, radii)

load M2;     %��ȡģ��,7*7��
%�������

Z00=conv2(M00,img);
Z11I=conv2(M11I,img);
Z11R=conv2(M11R,img);
Z20=conv2(M20,img);
Z31I=conv2(M31I,img);
Z31R=conv2(M31R,img);
Z40=conv2(M40,img);
%�ص����ಿ��
Z00=Z00(4:end-3,4:end-3);
Z11I=Z11I(4:end-3,4:end-3);
Z11R=Z11R(4:end-3,4:end-3);
Z20=Z20(4:end-3,4:end-3);
Z31I=Z31I(4:end-3,4:end-3);
Z31R=Z31R(4:end-3,4:end-3);
Z40=Z40(4:end-3,4:end-3);
%����k��oimg��theta��Zz11��l�ȵĳ�ʼ�����Է��������������С��ƥ�䡣
k=zeros(size(img,1),size(img,2));
oimg=k;
theta=k;
theta31=k;
Zz11=k;
Zz31=k;
z11=k;
z31=k;
Zz40=k;
Zz20=k;
h=k;
bounder=k;
l=100*ones(size(img,1),size(img,2));
l1=100*ones(size(img,1),size(img,2));
l2=100*ones(size(img,1),size(img,2));
%         ��ֹ���ӡ���ĸ��Ϊ0��������NAN
a=abs(Z11R)<0.0001;
b=abs(Z11I)<0.0001;
c=a&b;Z11R(c)=1;Z11I(c)=1;
theta=atan(Z11I./Z11R);
g=abs(Z31R)<0.0001;
j=abs(Z31I)<0.0001;
f=g&j;Z31R(f)=1;Z31I(f)=1;
theta31=atan(Z31I./Z31R);

Zz11=Z11R.*cos(theta)+Z11I.*sin(theta);
Zz31=Z31R.*cos(theta31)+Z31I.*sin(theta31); 
Zz20=Z20;
Zz40=Z40;
l1=sqrt((5*Zz40+3*Zz20)./(8*Zz20));
l2=sqrt((5*Zz31+Zz11)./(6*Zz11));
l=(l2+l1)/2;
k=1.5.*Zz11./((1-l2.^2).^1.5);
N=size(img,1);
[X,Y]=meshgrid(1:N,1:N);
x=X-(N+1)/2;y=(N+1)/2-Y;
index=(abs(l)<0.5)&(abs(k)>0.6);
xx=x+(7/2).*abs(l).*cos(theta31);
yy=y+(7/2).*abs(l).*sin(theta31);
for i=1:N
    for j=1:N
%         col(i,j)=round(xx(i,j)+(N+1)/2);
%         row(i,j)=round((N+1)/2-yy(i,j));
        col(i,j)=xx(i,j)+(N+1)/2;
        row(i,j)=round(N+1)/2-yy(i,j);
    end
end
imshow(img)
hold on
data=[];
for i=1:size(centers,1)
    point=[];
    for m=floor(centers(i,1)-radii(i))-1:ceil(centers(i,1)+radii(i))+1
        for n=floor(centers(i,2)-radii(i))-1:ceil(centers(i,2)+radii(i))+1
            if abs((sqrt((m-centers(i,1))^2+(n-centers(i,2))^2)-radii(i)))<0.707
                point=[point;[m,n,double(img(n,m))]];
            end
        end
    end
    data(i).point=point;
    plot(point(:,1),point(:,2),'r*');
end


end