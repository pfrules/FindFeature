% PixelPoints:检测钢球投影图像中心坐标，和对应帧数
function FirstOutput=FirstCalc(PixelPoints,pixel_scale)
[PsW, T_I_W,T_W_i, Pixel, Rd,R_i_I,Psi,PsI,init_angle] = initCalibrationGeometry(PixelPoints(1),pixel_scale,0,0);
gantry_angle=init_angle;
FirstOutput = [];
for i=1:size(PixelPoints,1)
%for i=1:2

    [PsW, T_I_W, T_W_i, Pixel, Rd,R_i_I,Psi,PsI,gantry_angle] = initCalibrationGeometry(PixelPoints(i),pixel_scale,gantry_angle,init_angle); 
    data.num=PixelPoints(i).num;
    data.PsW = PsW;
    data.T_I_W = T_I_W;
    data.Pixel = Pixel;
    data.Rd = Rd;
    data.theta= Rd(1);
    data.phi= Rd(2);
    data.eta= Rd(3);
    data.T_W_i=T_W_i;
    data.R_i_I=R_i_I;
    data.Psi=Psi;
    data.PsI=PsI; 
    data.init_angle=init_angle; 
    data.t=gantry_angle-init_angle;           
    FirstOutput = [FirstOutput; data];
	clc;
    disp(i);
end

end

%%
function [PsW, T_I_W,T_W_i, Pixel, Rd,R_i_I,Psi,PsI,gantry_angle] = initCalibrationGeometry(PixelPoints,pixel_scale,gantry_angle,init_angle)
%根据Ball Bearing标定模体，初步计算平板和X射线光源的位置和姿态，世界坐标系W与BB模体的自身坐标系重合。
%   输入： PixelPoints，        为BB模体中所有钢球球心在平板探测器上像素坐标，钢球需按顺序编号。
%                                  1号钢球初始位置位于世界坐标系W的X轴的正半轴,以及扫描帧数的编号的结构体。
%          init_angle，              体模初始角度估计,角度制
%          gantry_angle，         体模角度估计,角度制
%          pixel_scale,            物理尺寸像素尺寸比值
%
%   输出： T_I_W，                 4x4矩阵，现实平板坐标系到世界坐标系旋转矩阵
%          PsW，                 4x1向量，X光源相对于世界坐标系的位置向量，齐次坐标
%          Pixel,               2x1向量，世界坐标系原点在平板探测器上的投影点的像素坐标
%          Rd,                  3x1向量，对应于Cho论文中的phi，theta，eta，YXZ旋转顺序


%% 输入
XY = PixelPoints.point;
%输入检测
if(isempty(XY))
    warning('加载钢球球心像素文件错误！');
    PsW = -1;
    return;
end


%%
%求等中心投影点
%分组如下
% A组：（1-7，4-10）
% B组：（2-8，5-11）
% C组：（3-9，6-12）
w=[];
for i=1:6
    A = XY(i,1:2);
    B = XY(i+6,3:4);
    C = XY(i+6,1:2);
    D = XY(i,3:4); 
    c=GetIntersectPointofLines(A,B,C,D);
    if c==-1
        c=[];
    end
    w = [w;c];
end
center = [median(w(:,1)), median(w(:,2))];
Pixel = center';
 
% if (verbose)
%     figure;
%     plot(w(:,1),w(:,2), 'o'); %显示交点
%     hold on;
%     plot(center(:,1),center(:,2), '+k');
%     axis equal;
%     title('世界坐标系原点的投影点');
%     hold off;
% end

%%
%将所有椭圆的拟合点和穿透点转换到实际平板坐标系I
xyz1 = PixelToRealDetectorCoordinate(XY(:,1:2), center, pixel_scale);
xyz2 = PixelToRealDetectorCoordinate(XY(:,3:4), center, pixel_scale);

% if (verbose)
%     figure;
%     title('椭圆-平板坐标系I');
%     axis equal;
%     plot(xyz1(:,1),xyz1(:, 2),'+','Color','r');
%     hold on;
%     plot(xyz2(:,1),xyz2(:, 2),'+','Color','g');
%     hold off;
% end

%%
%拟合椭圆
[re1, cen1, ver1] = MyEllipseDirectFit(xyz1(:,1:2));
[re2, cen2, ver2] = MyEllipseDirectFit(xyz2(:,1:2));
[aa1,bb1]=getab(re1);
[aa2,bb2]=getab(re2);
% if (verbose)
%     line([ver1(1,1), ver1(2,1)], [ver1(1,2), ver1(2,2)], 'Color','r','LineStyle',':'); %绘制长轴1
%     line([ver2(1,1), ver2(2,1)], [ver2(1,2), ver2(2,2)], 'Color', 'g','LineStyle',':'); %绘制长轴2
% end

%%
%从平板理想坐标系i --> 平板实际坐标系I的旋转顺序为 YXZ，对应角度分别为 phi，theta，eta
%求解平板探测器的InPlaneAngle旋转角度eta
%重写椭圆1参数
%to rewrite the ellipse as p0*x^2+y^2-2*p1*x-2*p2*y+2*p3*x*y+p4=0
re1 = re1/re1(3);
p0=re1(1); p1=-re1(4)/2; p2=-re1(5)/2; p3=re1(2)/2; p4=re1(6);
%u=(p1-p2*p3)/(p0-p3^2);
%v=(p0*p2-p1*p3)/(p0-p3^2);
u1 = cen1(1,1); v1= cen1(1,2);
a1=p0/(p0*u1^2+v1^2+2*p3*u1*v1-p4);
b1=a1/p0;
c1=p3*b1;

%重写椭圆2参数
re2 = re2/re2(3);
p0=re2(1); p1=-re2(4)/2; p2=-re2(5)/2; p3=re2(2)/2; p4=re2(6);
%u=(p1-p2*p3)/(p0-p3^2);
%v=(p0*p2-p1*p3)/(p0-p3^2);
u2 = cen2(1,1); v2= cen2(1,2);
a2=p0/(p0*u2^2+v2^2+2*p3*u2*v2-p4);
b2=a2/p0;
c2=p3*b2;

p_0_1 = [u1, v1];
p_0_2 = [u2, v2];
 %此处使用论文中公式，转换坐标系后的结果将更准确，ab分别为椭圆参数
 p_0_a = (p_0_1*sqrt(a1/b1)+p_0_2*sqrt(a2/b2))/(sqrt(a1/b1)+sqrt(a2/b2));

%p_0_a = (p_0_1*sqrt(a2/b2)+p_0_2*sqrt(a1/b1))/(sqrt(a1/b1)+sqrt(a2/b2));
lp_a_1 = norm(p_0_1 - p_0_a);
lp_2_a = norm(p_0_a - p_0_2);
lp_1_2 = norm(p_0_1 - p_0_2);

x1 = ver1(1,1);x2 = ver1(2,1);y1 = ver1(1,2);y2 = ver1(2,2);
L1 = [x1-x2,y1-y2,0];
x1 = ver2(1,1);x2 = ver2(2,1);y1 = ver2(1,2);y2 = ver2(2,2);
L2 = [x1-x2,y1-y2,0];

angle1 = vectorAngleToX(L1);
angle2 = vectorAngleToX(L2);
eta = (lp_a_1 *angle1+lp_2_a*angle2)/lp_1_2; 
%%
%消除eta带来的影响，将数据点绕I坐标系原点旋转 -eta
r_xyz1 = RotateAlongZ(xyz1, -eta);
r_xyz2 = RotateAlongZ(xyz2, -eta);
%重新拟合椭圆
[r_re1, r_cen1, r_ver1] = MyEllipseDirectFit(r_xyz1(:,1:2));
[r_re2, r_cen2, r_ver2] = MyEllipseDirectFit(r_xyz2(:,1:2));
indxmin=find(r_xyz1(:,1)==min(r_xyz1(:,1)));
indxmax=find(r_xyz1(:,1)==max(r_xyz1(:,1)));
p_theta=GetIntersectPointofLines(r_xyz1(indxmin,1:2), r_xyz2(indxmin,1:2),r_xyz1(indxmax,1:2), r_xyz2(indxmax,1:2));
% if (verbose)
%     % 绘制椭圆长轴
%     hold on;
%     A=r_re1;
%     a=A(1);
%     b=A(2);
%     c=A(3);
%     d=A(4);
%     e=A(5);
%     f=A(6);
%     %eq0= 'a*x^2 + b*x*y + c*y^2 +d*x + e*y + f ';
%     eq0=@(x,y) a*x^2 + b*x*y + c*y^2 +d*x + e*y + f;
%eq0=@(x,y) a1*(x-r_cen1(1))^2 + b1*(y-cen1(2))^2 + 2*c1*(x-r_cen1(1))*(y-r_cen1(2))-1;
%     hold on;
%     h=ezplot(eq0,[-160,160,-180,180]);
%     set(h,'Color','k');
%     hold on;
%     line([r_ver1(1,1), r_ver1(2,1)], [r_ver1(1,2), r_ver1(2,2)], 'Color','r'); %绘制长轴1
%     line([r_ver2(1,1), r_ver2(2,1)], [r_ver2(1,2), r_ver2(2,2)], 'Color', 'g'); %绘制长轴2
% end

%重写椭圆1参数
re1 = r_re1/r_re1(3);
p0=re1(1); p1=-re1(4)/2; p2=-re1(5)/2; p3=re1(2)/2; p4=re1(6);
u1 = r_cen1(1,1); v1= r_cen1(1,2);
a1=p0/(p0*u1^2+v1^2+2*p3*u1*v1-p4);
b1=a1/p0;
c1=p3*b1;

%重写椭圆2参数
re2 = r_re2/r_re2(3);
p0=re2(1); p1=-re2(4)/2; p2=-re2(5)/2; p3=re2(2)/2; p4=re2(6);
u2 = r_cen2(1,1); v2= r_cen2(1,2);
a2=p0/(p0*u2^2+v2^2+2*p3*u2*v2-p4);
b2=a2/p0;
c2=p3*b2;

%%
%重新计算p0a
 p_0_1 = [u1, v1];
 p_0_2 = [u2, v2];
 %%此处使用自己新的公式，aa，bb是椭圆长短轴
 p_0_a = (p_0_1*(aa2/bb2)+p_0_2*(aa1/bb1))/((aa1/bb1)+(aa2/bb2));

%%
%求解theta，phi
%重新计算旋转之后的p_phi,标记为r_p_phi,重写椭圆拟合方程表达式
%首先求解Zsi
%采用椭圆的短轴近似的模拟L1和L2
radius = 115; %mm,钢球所组成圆的半径115
len = 210; %mm，两个钢球所组成的圆面之间的距离 210mm
L1=norm(r_ver1(3,:)-r_ver2(4,:));
L2=norm(r_ver1(4,:)-r_ver2(3,:));
Zsi = 2*radius*L1*L2/(len*(L2-L1));
Zdi=Zsi*(L1-len)/L1+radius;
r_p_phi=GetIntersectPointofLines(r_ver1(1,:), r_ver1(2,:),r_ver2(1,:), r_ver2(2,:));
if r_p_phi==-1
    phi=0;
else    
    eq=@(phi) sin(phi)+(c1*((r_p_phi(1)*sin(phi)*cos(phi))*a1*sqrt(a1)/...
    (sqrt(a1*b1+a1^2*b1*(r_p_phi(1)*sin(phi)*cos(phi))^2-c1^2)))/(2*a1)...
    -c2*((r_p_phi(1)*sin(phi)*cos(phi))*a2*sqrt(a2)/...
    (sqrt(a2*b2+a2^2*b2*(r_p_phi(1)*sin(phi)*cos(phi))^2-c2^2)))/(2*a2));
    options = optimset('Display','off','TolFun',1e-10,'TolX',1e-10);
    phi = fsolve(eq,0.01,options);
    phi=-sign(r_p_phi(1))*phi;
end
if p_theta==-1
theta=0;
else
theta=asin(Zsi*cos(phi)/p_theta(2));
end
%%

Rd = rad2deg([ theta,phi, -eta]'); %返回phi，theta，eta角度，弧度制

%%
%求解X光源位置，(X_s,Y_s,Z_s), 参考坐标系为平板实际坐标系I
Xsi=0; %由Cho等人定义的理想平板坐标系性质得出
Ysi = p_0_a(2)*cos(theta);

% R_I_i=(rotz(rad2deg(Rd(3)))*roty(rad2deg(Rd(2)))*rotx(rad2deg(Rd(1))));
R_I_i=(roty(Rd(2))*rotx(Rd(1))*rotz(Rd(3)));
R_i_I = R_I_i';





%%
%求解gantry的角度
len=105;
xyz = [ r_xyz1;r_xyz2]; %在实际平板坐标系I的XOY面内中对应于每个钢球的坐标
stepAngle = 30;
step = 0:1:11;
ang = deg2rad(step' * stepAngle);
xw = radius*sin(ang);
yw = -radius*cos(ang);
zw = len*ones(size(xw));
xyzw = [[xw, yw, zw];[xw, yw, -zw]];  % 在世界坐标系W中的钢球位置
xyzi=((roty(Rd(2))*rotx(Rd(1)))*xyz')';
[gantry_angle,Zsi,Zdi]=optimization(xyzi,xyzw,Zsi,Zdi,Ysi,gantry_angle);
%%
Psi = [Xsi, Ysi, Zsi]';
PsI = R_i_I * Psi;
%求解平板探测器的位置，(Ydi,Zdi)，参考坐标系为平板理想坐标系（Cho论文中的）
%根据Ford论文II.C.4小节的描述，推断(Y_d,Z_d)表示的世界坐标系原点在i坐标系中的坐标
%故，Xdi = Xwi，Ydi = Ywi， Zdi = Zwi
Xwi = 0;
Zwi = Zdi;
Ywi = Ysi*Zdi/Zsi;
Pwi = [Xwi, Ywi, Zwi]';
%%
%转换平板探测器的位置和方向矩阵、X光源位置到世界坐标系中
%转换平板
%R_i_W  = rotaY(-90)*rotaX(90)*rotaZ(t); % t为载物平台转角，与CBCT测试平台的转角方向相反
%t = rad2deg(angle); %载物平台转动角度 ,角度制, TODO:直接计算得到的角度不准确。
% t=rad2deg(gantry_angle-init_angle);
t=gantry_angle-init_angle;
 R_W_i = roty(-90)*rotx(90)*rotz(t);%世界坐标系到理想平板坐标系旋转矩阵
%R_W_i = rotz(90)*roty(-90);%世界坐标系到理想平板坐标系旋转矩阵
T_W_i = [R_W_i, Pwi; 0, 0, 0,1];%世界坐标系到理想平板坐标系变换矩阵
P_i_I = [0, 0, 0]';
T_i_I = [R_i_I, P_i_I; 0, 0, 0, 1]; %从平板理想坐标系i到平板实际坐标系I的变换矩阵
PsW = inv(T_W_i) * [Psi;1]; %X光源在世界坐标系W中的位置
% PsW=[roty(-t)*PsW(1:3);1];
T_I_W=inv(T_i_I*T_W_i);%现实平板坐标系到世界坐标系旋转矩阵

end

function [aa,bb]=getab(re)
a=re(1);
b=re(2)/2;
c=re(3);
d=re(4)/2;
f=re(5)/2;
g=re(6);
aa=sqrt((2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g))/((b^2-a*c)*(sqrt((a-c)^2+4*b^2)-(a+c))));
bb=sqrt((2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g))/((b^2-a*c)*(-sqrt((a-c)^2+4*b^2)-(a+c))));
if abs(aa)>abs(bb)
    temp=aa;
    aa=bb;
    bb=temp;
end
end


function P = GetIntersectPointofLines(A, B, C, D)
% GetIntersectPointofLines
%求平面内两直线AB和CD的交点
% http://www.cnblogs.com/DHUtoBUAA/

x1 = A(1);
y1 = A(2);
x2 = B(1);
y2 = B(2);
x3 = C(1);
y3 = C(2);
x4 = D(1);
y4 = D(2);

[A1, B1, C1] = GeneralEquation(x1,y1,x2,y2);
[A2, B2, C2] = GeneralEquation(x3,y3,x4,y4);

m=A1*B2-A2*B1;
if m==0
        disp('无交点！');
        P = -1; 
else
        x=(C2*B1-C1*B2)/m;
        y=(C1*A2-C2*A1)/m;
        P = [x, y];
end

end

function [A,B,C] = GeneralEquation(first_x,first_y,second_x,second_y)
    %一般式 Ax+By+C=0
    % from http://www.cnblogs.com/DHUtoBUAA/
    A=second_y-first_y;
    B=first_x-second_x;
    C=second_x*first_y-first_x*second_y;
end

function p = RotateAlongZ(point, angle)
%将像素坐标转换到实际平板坐标系中
%输入：   point, nx3,（n，1）x坐标，(n,2)y坐标，(n,3)z坐标
%         angle为绕Z轴的旋转角度,弧度制
%输出：   p 为nx3的矩阵，旋转后的点坐标
[m, n] = size(point);
t = point(:,1:2);

A = [cos(angle) -sin(angle);
    sin(angle), cos(angle)];

re = [];
for i=1:m
    s = A*(t(i,:))';
    re = [re; s'];
end

p = [re,point(:,3)];
end

function p = PixelToRealDetectorCoordinate(pixel, center, scale)
%将像素坐标转换到实际平板坐标系中
%输入：   pixel（i，1）平板的像素u轴坐标，pixel（i，2）位v轴坐标，i为行号
%        center(x,y)为穿透点像素坐标
%        scale为像素大小，单位为mm
%输出：   p 为nx3的矩阵，为pixel对应像素在实际平板坐标系中的坐标值
[m, n] = size(pixel);
w=[];
for i=1:m
    % 像素坐标系（u，v）与平板坐标系I（x,y,z）, u和x方向相反，v和y的方向相同
    t = [-(pixel(i,1)-center(1))/scale, (pixel(i,2)-center(2))/scale];
    w = [w;t];
end
z = zeros(m,1);
p=[w,z];
end


function  ang = vectorAngleToX(v)
%求解输入向量v与X轴的夹角大小与方向
%   输入： v，                  1x3，输入的向量。
%
%   输出： ang，                标量，弧度制，v与X轴的夹角，逆时针为正，取值范围[-pi/2,pi/2]。

La = v;
n_x = [1,0,0];
c = cross(n_x, La);
angle = acos(dot(La, n_x)/(norm(La)*norm(n_x)));

if ( c(3)>0 )
    if ((0< angle) && (angle< pi/2))
        eta = angle;
    else
        eta = angle-pi;
    end
elseif ( c(3)<0 )
    if ((0< angle) && (angle< pi/2))
        eta = - angle;
    else
        eta = pi-angle;
    end
else
    disp('请注意，包含零向量求角度！');
    eta = 0;
end

% n_z = [0,0,1];
% angle_d = acos(dot(c,n_z)/(norm(c)*norm(n_z))) % 判断角度方向，逆时针为正。
% if ((angle_d >=0 ) && (angle_d < pi/2)) % 正向
%     if ((0< angle) && (angle< pi/2))
%         eta = angle;
%     else
%         eta = angle-pi;
%     end
% else                                   %负向
%     if ((0< angle) && (angle< pi/2))
%         eta = - angle;
%     else
%         eta = pi-angle;
%     end
% end

ang = eta;

end


function [Re, center, vertex]= MyEllipseDirectFit(XY)
%  Direct ellipse fit, proposed in article
%    A. W. Fitzgibbon, M. Pilu, R. B. Fisher
%     "Direct Least Squares Fitting of Ellipses"
%     IEEE Trans. PAMI, Vol. 21, pages 476-480 (1999)
%
%  Our code is based on a numerically stable version
%  of this fit published by R. Halir and J. Flusser
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: A = [a b c d e f]' is the vector of algebraic 
%             parameters of the fitting ellipse:
%             ax^2 + bxy + cy^2 +dx + ey + f = 0
%             the vector A is normed, so that ||A||=1
%
%  This is a fast non-iterative ellipse fit.
%
%  It returns ellipses only, even if points are
%  better approximated by a hyperbola.
%  It is somewhat biased toward smaller ellipses.
%
centroid = mean(XY); % the centroid of the data set
warning off;
D1 = [(XY(:,1)-centroid(1)).^2, (XY(:,1)-centroid(1)).*(XY(:,2)-centroid(2)),...
      (XY(:,2)-centroid(2)).^2];
D2 = [XY(:,1)-centroid(1), XY(:,2)-centroid(2), ones(size(XY,1),1)];
S1 = D1'*D1;
S2 = D1'*D2;
S3 = D2'*D2;
T = -inv(S3)*S2';
M = S1 + S2*T;
M = [M(3,:)./2; -M(2,:); M(1,:)./2];
[evec,eval] = eig(M);
cond = 4*evec(1,:).*evec(3,:)-evec(2,:).^2;
A1 = evec(:,find(cond>0));
A = [A1; T*A1];
A4 = A(4)-2*A(1)*centroid(1)-A(2)*centroid(2);
A5 = A(5)-2*A(3)*centroid(2)-A(2)*centroid(1);
A6 = A(6)+A(1)*centroid(1)^2+A(3)*centroid(2)^2+...
     A(2)*centroid(1)*centroid(2)-A(4)*centroid(1)-A(5)*centroid(2);
A(4) = A4;  A(5) = A5;  A(6) = A6;
% if A(1)<0
%     A=-A;
% end
Re=A;

%长轴倾角
%ellipse equation: a*x^2 + b*x*y + c*y^2 +d*x + e*y + 1 = 0;
% Refer to：http://blog.csdn.net/ningyaliuhebei/article/details/46327681
B=Re(2)/Re(6); A=Re(1)/Re(6);C=Re(3)/Re(6); D=Re(4)/Re(6);E=Re(5)/Re(6);
if B==0
     if A<C/2
    theta1=0; 
    else
      theta1= pi/2;
    end  
else
    if A<C/2
    theta1=0.5*atan(B/(A-C)); 
    else
      theta1= pi/2 + 0.5*atan(B/(A-C));
    end
end
K1= tan(theta1); %长轴斜率
if (abs(K1)<1e-8) 
    K1=0;
end
%椭圆的中心点
xc=(B*E-2*C*D)/(4*A*C-B^2);
yc=(B*D-2*A*E)/(4*A*C-B^2);
%椭圆的长、短半轴长度
la=sqrt(2*(A*xc^2+C*yc^2+B*xc*yc-1)/(A+C-sqrt((A-C)^2+B^2)));
lb=sqrt(2*(A*xc^2+C*yc^2+B*xc*yc-1)/(A+C+sqrt((A-C)^2+B^2)));
%长、短轴顶点坐标
n1=[1,K1]; n1=n1/norm(n1);
center = [xc, yc];
va1=center+la*n1;  %长轴顶点1, 坐标轴正向
va2=center-la*n1;  %长轴顶点2，坐标轴负向

if(K1 == 0)
    n2=[0,1];
else
    n2=[1,-1/K1];
end
n2=n2/norm(n2);
if n2(2) >= 0
    vb1=center+lb*n2;  %短轴顶点1，坐标轴正向
    vb2=center-lb*n2;  %短轴顶点2，坐标轴负向
else
    vb1=center-lb*n2;  %短轴顶点1，坐标轴正向
    vb2=center+lb*n2;  %短轴顶点2，坐标轴负向
end

center = [xc, yc];
vertex = [va1;va2;vb1;vb2];

end  %  EllipseDirectFit



function [gantry_angle,Zsi,Zdi]=optimization(xyz,xyzw,Zsi,Zdi,Ysi,gantry_angle)
    para=[gantry_angle,Zsi,Zdi];
%     options = optimset('LargeScale','off','LevenbergMarquardt','on');
       options = optimset('Display','off','TolX',eps,'TolFun',eps,'LargeScale','off','Algorithm','Levenberg-Marquardt');
    [x,resnorm,residual,exitflag,output]  = lsqnonlin( @myop, para, [],[],options, xyz, xyzw,Ysi);
%  options = optimset('Display','off','TolX',eps,'TolFun',eps,'LargeScale','off','Algorithm','trust-region-reflective');
%   [x,resnorm,residual,exitflag,output]  = lsqnonlin( @myop, para, [-2*pi/180,-2*pi/180,-pi],[5*pi/180,5*pi/180,pi],options, xyz, xyzw, Zsi,Zdi,oa);
     % display the result
    gantry_angle=x(1);
    Zsi=x(2);
    Zdi=x(3);
end


function f = myop(para, xyz, xyzw, Ysi)
    t=para(1);
    Zsi=para(2);
    Zdi=para(3);
    xyzi=xyz+repmat([0,-Ysi,-Zsi],size(xyz,1),1);
    xyz0=(roty(-90)*rotx(90)*rotz(t)*xyzw')'+repmat([0,-Ysi*(Zsi-Zdi)/Zsi,Zdi-Zsi],size(xyz,1),1);
    a=xyzi(:,3).*xyz0(:,2)-xyzi(:,2).*xyz0(:,3);
    b=xyzi(:,1).*xyz0(:,3)-xyzi(:,3).*xyz0(:,1);
%   c=xyzi(:,2).*xyz0(:,1)-xyzi(:,1).*xyz0(:,2);
 f=[a,b];
%f=[a,b];
% f=a+b;
% f=a.^2+b.^2+c.^2;
end