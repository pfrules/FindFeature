function result=lowway(img,centers, radii)
%% zernike方法
% load M2;     %获取模板,7*7的
%卷积运算

% Z00=conv2(M00,img);
% Z11I=conv2(M11I,img);
% Z11R=conv2(M11R,img);
% Z20=conv2(M20,img);
% Z31I=conv2(M31I,img);
% Z31R=conv2(M31R,img);
% Z40=conv2(M40,img);
% %截掉多余部分
% Z00=Z00(4:end-3,4:end-3);
% Z11I=Z11I(4:end-3,4:end-3);
% Z11R=Z11R(4:end-3,4:end-3);
% Z20=Z20(4:end-3,4:end-3);
% Z31I=Z31I(4:end-3,4:end-3);
% Z31R=Z31R(4:end-3,4:end-3);
% Z40=Z40(4:end-3,4:end-3);
% %设置k，oimg，theta，Zz11，l等的初始矩阵，以防后来矩阵运算大小不匹配。
% k=zeros(size(img,1),size(img,2));
% oimg=k;
% theta=k;
% theta31=k;
% Zz11=k;
% Zz31=k;
% z11=k;
% z31=k;
% Zz40=k;
% Zz20=k;
% h=k;
% bounder=k;
% l=100*ones(size(img,1),size(img,2));
% l1=100*ones(size(img,1),size(img,2));
% l2=100*ones(size(img,1),size(img,2));
% %         防止分子、分母均为0而产生的NAN
% a=abs(Z11R)<0.0001;
% b=abs(Z11I)<0.0001;
% c=a&b;Z11R(c)=1;Z11I(c)=1;
% theta=atan(Z11I./Z11R);
% g=abs(Z31R)<0.0001;
% j=abs(Z31I)<0.0001;
% f=g&j;Z31R(f)=1;Z31I(f)=1;
% theta31=atan(Z31I./Z31R);
% 
% Zz11=Z11R.*cos(theta)+Z11I.*sin(theta);
% Zz31=Z31R.*cos(theta31)+Z31I.*sin(theta31); 
% Zz20=Z20;
% Zz40=Z40;
% l1=sqrt((5*Zz40+3*Zz20)./(8*Zz20));
% l2=sqrt((5*Zz31+Zz11)./(6*Zz11));
% l=(l2+l1)/2;
% k=1.5.*Zz11./((1-l2.^2).^1.5);
% N=size(img,1);
% [X,Y]=meshgrid(1:N,1:N);
% x=X-(N+1)/2;y=(N+1)/2-Y;
% index=(abs(l)<0.5)&(abs(k)>0.6);
% xx=x+(7/2).*abs(l).*cos(theta31);
% yy=y+(7/2).*abs(l).*sin(theta31);
% for i=1:N
%     for j=1:N
% %         col(i,j)=round(xx(i,j)+(N+1)/2);
% %         row(i,j)=round((N+1)/2-yy(i,j));
%         col(i,j)=xx(i,j)+(N+1)/2;
%         row(i,j)=round(N+1)/2-yy(i,j);
%     end
% end

%%
%线性拟合
% imshow(img)
% hold on
data=[];
for i=1:size(centers,1)
    point=[];
    for m=floor(centers(i,1)-radii(i))-1:ceil(centers(i,1)+radii(i))+1
        for n=floor(centers(i,2)-radii(i))-1:ceil(centers(i,2)+radii(i))+1
            dis=sqrt((m-centers(i,1))^2+(n-centers(i,2))^2);
            if abs((dis-radii(i)))<0.707
                point=[point;[m,n,double(img(n,m)),dis]];
            end
        end
    end
    data(i).point=point;
   % plot(point(:,1),point(:,2),'r*');
end
for p=1:size(centers,1)
point=data(p).point;
% [mu,sigma]=normfit(point(:,3));
% [y,x]=hist(point(:,3),20);
% M=[];
% M(:,1)=x.^2;
% M(:,2)=x;
% M(:,3)=1;
% M(:,4)=-y';
% [U,D,V]=svd(M);
% V=V(1:end-1,end)./V(end,end);
% mu=-V(2)/(2*V(1));
% % if mu>max(point(:,3))||mu<min(point(:,3))
% %     mu=max(point(:,3));
% % end
% mu=x(y==max(y));
% mu=mu(end);
% %% 画直方图
% 
% 
% figure
% bar(x,y,'FaceColor','r','EdgeColor','w');box off
% xlim([mu-3*sigma,mu+3*sigma])
% a2=axes;
% % %ezplot(@(x)normpdf(x,mu,sigma),[mu-3*sigma,mu+3*sigma])
% % ezplot(@(grayscale)V(1)*grayscale*grayscale+V(2)*grayscale+V(3),[min(x),max(x)])
% set(a2,'box','off','yaxislocation','right','color','none')
% title '频数直方图与正态分布密度函数（拟合）'
%%
m=size(point,1);                     % 点数
A={zeros(0,2)};                  % 元包数组中仅包含一个元素
A_g_v=repmat(A,m,1);             % 保存灰度值
A_p=repmat(A,m,1);               % 保存5*5坐标值
Coef=repmat(A,m,1);              % 系数
Newpoint=zeros(m,2);
D=repmat(A,25,1);  
for i=1:m
    x=point(i,1);
    y=point(i,2);
    A_g_v{i}=zeros(5,5);
    kk=1;
    for j=-2:2
        for k=-2:2
            A_p{i}(kk,1)=j;                % X方向
            A_p{i}(kk,2)=k;                % Y方向            
            A_p{i}(kk,3)=img(y+k,x+j);      % 灰度值 向量表示
            A_g_v{i}(k+3,j+3)=img(y+k,x+j); % 灰度值 矩阵表示          
            kk=kk+1;
        end
    end
end
%[XX,YY]=meshgrid(A_p{point}(1,1):1:A_p{point}(25,1),...
%                 A_p{point}(1,2):1:A_p{point}(25,2));
%ZZ=A_g_v{point};
%mesh(XX,YY,ZZ);

for j=1:25
    D{j}=zeros(10,1);
    D{j}(1)=1;
    D{j}(2)=A_p{1}(j,1);
    D{j}(3)=A_p{1}(j,2);
    D{j}(4)=A_p{1}(j,1)^2;
    D{j}(5)=A_p{1}(j,1)*A_p{1}(j,2);
    D{j}(6)=A_p{1}(j,2)^2;
    D{j}(7)=A_p{1}(j,1)^3;
    D{j}(8)=A_p{1}(j,1)^2*A_p{1}(j,2);
    D{j}(9)=A_p{1}(j,1)*A_p{1}(j,2)^2;
    D{j}(10)=A_p{1}(j,2)^3;
end
K=[D{1}(1) D{1}(2) D{1}(3) D{1}(4) D{1}(5) D{1}(6) D{1}(7) D{1}(8) D{1}(9) D{1}(10);
   D{2}(1) D{2}(2) D{2}(3) D{2}(4) D{2}(5) D{2}(6) D{2}(7) D{2}(8) D{2}(9) D{2}(10);
   D{3}(1) D{3}(2) D{3}(3) D{3}(4) D{3}(5) D{3}(6) D{3}(7) D{3}(8) D{3}(9) D{3}(10);
   D{4}(1) D{4}(2) D{4}(3) D{4}(4) D{4}(5) D{4}(6) D{4}(7) D{4}(8) D{4}(9) D{4}(10);
   D{5}(1) D{5}(2) D{5}(3) D{5}(4) D{5}(5) D{5}(6) D{5}(7) D{5}(8) D{5}(9) D{5}(10);
   D{6}(1) D{6}(2) D{6}(3) D{6}(4) D{6}(5) D{6}(6) D{6}(7) D{6}(8) D{6}(9) D{6}(10);
   D{7}(1) D{7}(2) D{7}(3) D{7}(4) D{7}(5) D{7}(6) D{7}(7) D{7}(8) D{7}(9) D{7}(10);
   D{8}(1) D{8}(2) D{8}(3) D{8}(4) D{8}(5) D{8}(6) D{8}(7) D{8}(8) D{8}(9) D{8}(10);
   D{9}(1) D{9}(2) D{9}(3) D{9}(4) D{9}(5) D{9}(6) D{9}(7) D{9}(8) D{9}(9) D{9}(10);
   D{10}(1) D{10}(2) D{10}(3) D{10}(4) D{10}(5) D{10}(6) D{10}(7) D{10}(8) D{10}(9) D{10}(10);
   D{11}(1) D{11}(2) D{11}(3) D{11}(4) D{11}(5) D{11}(6) D{11}(7) D{11}(8) D{11}(9) D{11}(10);
   D{12}(1) D{12}(2) D{12}(3) D{12}(4) D{12}(5) D{12}(6) D{12}(7) D{12}(8) D{12}(9) D{12}(10);
   D{13}(1) D{13}(2) D{13}(3) D{13}(4) D{13}(5) D{13}(6) D{13}(7) D{13}(8) D{13}(9) D{13}(10);
   D{14}(1) D{14}(2) D{14}(3) D{14}(4) D{14}(5) D{14}(6) D{14}(7) D{14}(8) D{14}(9) D{14}(10);
   D{15}(1) D{15}(2) D{15}(3) D{15}(4) D{15}(5) D{15}(6) D{15}(7) D{15}(8) D{15}(9) D{15}(10);
   D{16}(1) D{16}(2) D{16}(3) D{16}(4) D{16}(5) D{16}(6) D{16}(7) D{16}(8) D{16}(9) D{16}(10);
   D{17}(1) D{17}(2) D{17}(3) D{17}(4) D{17}(5) D{17}(6) D{17}(7) D{17}(8) D{17}(9) D{17}(10);
   D{18}(1) D{18}(2) D{18}(3) D{18}(4) D{18}(5) D{18}(6) D{18}(7) D{18}(8) D{18}(9) D{18}(10);
   D{19}(1) D{19}(2) D{19}(3) D{19}(4) D{19}(5) D{19}(6) D{19}(7) D{19}(8) D{19}(9) D{19}(10);
   D{20}(1) D{20}(2) D{20}(3) D{20}(4) D{20}(5) D{20}(6) D{20}(7) D{20}(8) D{20}(9) D{20}(10);
   D{21}(1) D{21}(2) D{21}(3) D{21}(4) D{21}(5) D{21}(6) D{21}(7) D{21}(8) D{21}(9) D{21}(10);
   D{22}(1) D{22}(2) D{22}(3) D{22}(4) D{22}(5) D{22}(6) D{22}(7) D{22}(8) D{22}(9) D{22}(10);
   D{23}(1) D{23}(2) D{23}(3) D{23}(4) D{23}(5) D{23}(6) D{23}(7) D{23}(8) D{23}(9) D{23}(10);
   D{24}(1) D{24}(2) D{24}(3) D{24}(4) D{24}(5) D{24}(6) D{24}(7) D{24}(8) D{24}(9) D{24}(10);
   D{25}(1) D{25}(2) D{25}(3) D{25}(4) D{25}(5) D{25}(6) D{25}(7) D{25}(8) D{25}(9) D{25}(10);];

sin_Theta=zeros(m,1);
cos_Theta=zeros(m,1);
r=zeros(m,1);
%T=zeros(m,4);

% 梯度最大化
for i=1:m   
    KK=Coefficient(A_p,i);
    Coef{i}=zeros(10,1);
    Coef{i}=K\KK;                   % det(D)表示求D的行列式,cramer法则求解
    K2=Coef{i}(2,1);                % 
    K3=Coef{i}(3,1);                % 
    cos_Theta(i)=K2/sqrt(K2^2+K3^2);%Theta表示梯度方向与X轴夹角
    sin_Theta(i)=K3/sqrt(K2^2+K3^2);

    a=6*(Coef{i}(7,1)*(sin_Theta(i)^3)+Coef{i}(8,1)*(sin_Theta(i)^2)*cos_Theta(i)+...
        Coef{i}(9,1)*sin_Theta(i)*(cos_Theta(i)^2)+Coef{i}(10,1)*(cos_Theta(i)^3));
%     a=6*(Coef{i}(7,1)*(sin_Theta(i)^3)+Coef{i}(10,1)*(cos_Theta(i)^3))+...
%     2*(Coef{i}(8,1)*sin_Theta(i)*(1+cos_Theta(i)^2)+...
%     Coef{i}(9,1)*cos_Theta(i)*(1+sin_Theta(i)^2));
    b=2*(Coef{i}(4,1)*(sin_Theta(i)^2)+Coef{i}(5,1)*sin_Theta(i)*cos_Theta(i)+...
        Coef{i}(6,1)*(cos_Theta(i)^2));
    
    r(i)=-b/a; %出问题了？
    if sin_Theta(i)*cos_Theta(i)<0
        r(i)=b/a; 
    end
    
%     if abs(cos_Theta(i))>0.99
%          T(i,4)=0;
%     elseif (r(i)>0&&sin_Theta(i)>0&&cos_Theta(i)>0)||(r(i)<0&&sin_Theta(i)<0&&cos_Theta(i)<0)
%         T(i,4)=1;
%     elseif (r(i)>0&&sin_Theta(i)>0&&cos_Theta(i)<0)||(r(i)<0&&sin_Theta(i)<0&&cos_Theta(i)>0)
%         T(i,4)=2;
%     elseif (r(i)>0&&sin_Theta(i)<0&&cos_Theta(i)>0)||(r(i)<0&&sin_Theta(i)>0&&cos_Theta(i)<0)
%         T(i,4)=4;
%     else
%         T(i,4)=3;    
%     end
%     T(i,1)=r(i);T(i,2)=sin_Theta(i);T(i,3)=cos_Theta(i);
    
    if abs(r(i))>=2
        r(i)=r(i-1);
    end

    Newpoint(i,1)=point(i,1)+r(i)*cos_Theta(i);      % x=x0+r*cos(Theta);
    Newpoint(i,2)=point(i,2)+r(i)*sin_Theta(i);      % y=y0+r*sin(Theta);
end


%% 灰度值相等法
% for i=1:m   
%     KK=Coefficient(A_p,i);
%     Coef{i}=zeros(10,1);
%     %原作者方法
%     %Coef{i}=K\KK;                   % det(D)表示求D的行列式,cramer法则求解
%     %我的方法
%     M=[K,-KK];
%     [U,D,V]=svd(M);
%     V=V/V(end,end);
%     Coef{i}=V(1:end-1,end);
%     
%     K2=Coef{i}(2,1);                % 
%     K3=Coef{i}(3,1);                % 
%     %原作者
% %     sin_Theta(i)=K2/sqrt(K2^2+K3^2);%Theta表示梯度方向与X轴夹角
% %     cos_Theta(i)=K3/sqrt(K2^2+K3^2);
% 
%     %修正
%     cos_Theta(i)=K2/sqrt(K2^2+K3^2);%Theta表示梯度方向与X轴夹角
%     sin_Theta(i)=K3/sqrt(K2^2+K3^2);
%     
%     a=Coef{i}(1,1)-mu;
%     b=Coef{i}(2,1)*cos_Theta(i)+Coef{i}(3,1)*sin_Theta(i);
%     c=Coef{i}(4,1)*cos_Theta(i)^2+Coef{i}(5,1)*cos_Theta(i)*sin_Theta(i)+Coef{i}(6,1)*sin_Theta(i)^2;
%     d=Coef{i}(7,1)*cos_Theta(i)^3+Coef{i}(8,1)*cos_Theta(i)^2*sin_Theta(i)+Coef{i}(9,1)*cos_Theta(i)*sin_Theta(i)^2+Coef{i}(10,1)*sin_Theta(i)^3;
%     options = optimset('Display','off','TolFun',1e-10,'TolX',1e-10);
%     eq=@(rol) a+b*rol+c*rol^2+d*rol^3;
%     rol = fsolve(eq,0.01,options);
%     %figure;ezplot(@(rol) a+b*rol+c*rol^2+d*rol^3,[-3,3])
%     r(i)=rol; 
%     Newpoint(i,1)=point(i,1)+r(i)*cos_Theta(i);      % x=x0+r*cos(Theta);
%     Newpoint(i,2)=point(i,2)+r(i)*sin_Theta(i);      % y=y0+r*sin(Theta);
%     kkk=Coef{i};
%     xx=r(i)*cos_Theta(i);
%     yy=r(i)*sin_Theta(i);
%     Newpoint(i,3)=kkk(1)+kkk(2)*xx+kkk(3)*yy+kkk(4)*xx*xx+kkk(5)*yy*xx+kkk(6)*yy*yy+...
%         kkk(7)*xx^3+kkk(8)*xx^2*yy+kkk(9)*xx*yy^2+kkk(10)*yy^3;
%     Newpoint(i,4)=a+b*rol+c*rol^2+d*rol^3;
% end
% endpoint=[];
% for i=1:m
%     if abs(Newpoint(i,3)-mu)>1&&abs(Newpoint(i,4))>1
%         continue;
%     end
%     endpoint=[endpoint;Newpoint(i,:)];
% end
%%
%figure(p);%imshow(img);
hold on
plot(Newpoint(:,1),Newpoint(:,2),'k*');                 %修正前

%figure;
%hold on
%plot(endpoint(:,1),endpoint(:,2),'b*');          %修正后
% plot(point(:,1),point(:,2),'r*');          %修正后

% result.point=endpoint;
[Re, center2, vertex]= MyEllipseDirectFit(Newpoint(:,1:2));
    a=Re(1);
    b=Re(2);
    c=Re(3);
    d=Re(4);
    e=Re(5);
    f=Re(6);
    %eq0= 'a*x^2 + b*x*y + c*y^2 +d*x + e*y + f ';
    eq0=@(x,y) a*x^2 + b*x*y + c*y^2 +d*x + e*y + f;
    hold on;
    h=ezplot(eq0,[center2(1)-14,center2(1)+14,center2(2)-14,center2(2)+14]);
    set(h,'Color','b');
% plot(center2(1),center2(2),'k*');
% plot(centers(p,1),centers(p,2),'r*');
% result.Re=Re;
% result.center=center;
% result.vertex=vertex;
% resultdata=[resultdata;result];
result(p,:)=center2;
end

end



function B=Coefficient(A_p,i)
B=[A_p{i}(1,3);A_p{i}(2,3);A_p{i}(3,3);A_p{i}(4,3);A_p{i}(5,3);
    A_p{i}(6,3);A_p{i}(7,3);A_p{i}(8,3);A_p{i}(9,3);A_p{i}(10,3);
    A_p{i}(11,3);A_p{i}(12,3);A_p{i}(13,3);A_p{i}(14,3);A_p{i}(15,3);
    A_p{i}(16,3);A_p{i}(17,3);A_p{i}(18,3);A_p{i}(19,3);A_p{i}(20,3);
    A_p{i}(21,3);A_p{i}(22,3);A_p{i}(23,3);A_p{i}(24,3);A_p{i}(25,3)];
end