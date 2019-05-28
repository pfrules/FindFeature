% PixelPoints:������ͶӰͼ���������꣬�Ͷ�Ӧ֡��
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
%����Ball Bearing�궨ģ�壬��������ƽ���X���߹�Դ��λ�ú���̬����������ϵW��BBģ�����������ϵ�غϡ�
%   ���룺 PixelPoints��        ΪBBģ�������и���������ƽ��̽�������������꣬�����谴˳���š�
%                                  1�Ÿ����ʼλ��λ����������ϵW��X���������,�Լ�ɨ��֡���ı�ŵĽṹ�塣
%          init_angle��              ��ģ��ʼ�Ƕȹ���,�Ƕ���
%          gantry_angle��         ��ģ�Ƕȹ���,�Ƕ���
%          pixel_scale,            ����ߴ����سߴ��ֵ
%
%   ����� T_I_W��                 4x4������ʵƽ������ϵ����������ϵ��ת����
%          PsW��                 4x1������X��Դ�������������ϵ��λ���������������
%          Pixel,               2x1��������������ϵԭ����ƽ��̽�����ϵ�ͶӰ�����������
%          Rd,                  3x1��������Ӧ��Cho�����е�phi��theta��eta��YXZ��ת˳��


%% ����
XY = PixelPoints.point;
%������
if(isempty(XY))
    warning('���ظ������������ļ�����');
    PsW = -1;
    return;
end


%%
%�������ͶӰ��
%��������
% A�飺��1-7��4-10��
% B�飺��2-8��5-11��
% C�飺��3-9��6-12��
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
%     plot(w(:,1),w(:,2), 'o'); %��ʾ����
%     hold on;
%     plot(center(:,1),center(:,2), '+k');
%     axis equal;
%     title('��������ϵԭ���ͶӰ��');
%     hold off;
% end

%%
%��������Բ����ϵ�ʹ�͸��ת����ʵ��ƽ������ϵI
xyz1 = PixelToRealDetectorCoordinate(XY(:,1:2), center, pixel_scale);
xyz2 = PixelToRealDetectorCoordinate(XY(:,3:4), center, pixel_scale);

% if (verbose)
%     figure;
%     title('��Բ-ƽ������ϵI');
%     axis equal;
%     plot(xyz1(:,1),xyz1(:, 2),'+','Color','r');
%     hold on;
%     plot(xyz2(:,1),xyz2(:, 2),'+','Color','g');
%     hold off;
% end

%%
%�����Բ
[re1, cen1, ver1] = MyEllipseDirectFit(xyz1(:,1:2));
[re2, cen2, ver2] = MyEllipseDirectFit(xyz2(:,1:2));
[aa1,bb1]=getab(re1);
[aa2,bb2]=getab(re2);
% if (verbose)
%     line([ver1(1,1), ver1(2,1)], [ver1(1,2), ver1(2,2)], 'Color','r','LineStyle',':'); %���Ƴ���1
%     line([ver2(1,1), ver2(2,1)], [ver2(1,2), ver2(2,2)], 'Color', 'g','LineStyle',':'); %���Ƴ���2
% end

%%
%��ƽ����������ϵi --> ƽ��ʵ������ϵI����ת˳��Ϊ YXZ����Ӧ�Ƕȷֱ�Ϊ phi��theta��eta
%���ƽ��̽������InPlaneAngle��ת�Ƕ�eta
%��д��Բ1����
%to rewrite the ellipse as p0*x^2+y^2-2*p1*x-2*p2*y+2*p3*x*y+p4=0
re1 = re1/re1(3);
p0=re1(1); p1=-re1(4)/2; p2=-re1(5)/2; p3=re1(2)/2; p4=re1(6);
%u=(p1-p2*p3)/(p0-p3^2);
%v=(p0*p2-p1*p3)/(p0-p3^2);
u1 = cen1(1,1); v1= cen1(1,2);
a1=p0/(p0*u1^2+v1^2+2*p3*u1*v1-p4);
b1=a1/p0;
c1=p3*b1;

%��д��Բ2����
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
 %�˴�ʹ�������й�ʽ��ת������ϵ��Ľ������׼ȷ��ab�ֱ�Ϊ��Բ����
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
%����eta������Ӱ�죬�����ݵ���I����ϵԭ����ת -eta
r_xyz1 = RotateAlongZ(xyz1, -eta);
r_xyz2 = RotateAlongZ(xyz2, -eta);
%���������Բ
[r_re1, r_cen1, r_ver1] = MyEllipseDirectFit(r_xyz1(:,1:2));
[r_re2, r_cen2, r_ver2] = MyEllipseDirectFit(r_xyz2(:,1:2));
indxmin=find(r_xyz1(:,1)==min(r_xyz1(:,1)));
indxmax=find(r_xyz1(:,1)==max(r_xyz1(:,1)));
p_theta=GetIntersectPointofLines(r_xyz1(indxmin,1:2), r_xyz2(indxmin,1:2),r_xyz1(indxmax,1:2), r_xyz2(indxmax,1:2));
% if (verbose)
%     % ������Բ����
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
%     line([r_ver1(1,1), r_ver1(2,1)], [r_ver1(1,2), r_ver1(2,2)], 'Color','r'); %���Ƴ���1
%     line([r_ver2(1,1), r_ver2(2,1)], [r_ver2(1,2), r_ver2(2,2)], 'Color', 'g'); %���Ƴ���2
% end

%��д��Բ1����
re1 = r_re1/r_re1(3);
p0=re1(1); p1=-re1(4)/2; p2=-re1(5)/2; p3=re1(2)/2; p4=re1(6);
u1 = r_cen1(1,1); v1= r_cen1(1,2);
a1=p0/(p0*u1^2+v1^2+2*p3*u1*v1-p4);
b1=a1/p0;
c1=p3*b1;

%��д��Բ2����
re2 = r_re2/r_re2(3);
p0=re2(1); p1=-re2(4)/2; p2=-re2(5)/2; p3=re2(2)/2; p4=re2(6);
u2 = r_cen2(1,1); v2= r_cen2(1,2);
a2=p0/(p0*u2^2+v2^2+2*p3*u2*v2-p4);
b2=a2/p0;
c2=p3*b2;

%%
%���¼���p0a
 p_0_1 = [u1, v1];
 p_0_2 = [u2, v2];
 %%�˴�ʹ���Լ��µĹ�ʽ��aa��bb����Բ������
 p_0_a = (p_0_1*(aa2/bb2)+p_0_2*(aa1/bb1))/((aa1/bb1)+(aa2/bb2));

%%
%���theta��phi
%���¼�����ת֮���p_phi,���Ϊr_p_phi,��д��Բ��Ϸ��̱��ʽ
%�������Zsi
%������Բ�Ķ�����Ƶ�ģ��L1��L2
radius = 115; %mm,���������Բ�İ뾶115
len = 210; %mm��������������ɵ�Բ��֮��ľ��� 210mm
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

Rd = rad2deg([ theta,phi, -eta]'); %����phi��theta��eta�Ƕȣ�������

%%
%���X��Դλ�ã�(X_s,Y_s,Z_s), �ο�����ϵΪƽ��ʵ������ϵI
Xsi=0; %��Cho���˶��������ƽ������ϵ���ʵó�
Ysi = p_0_a(2)*cos(theta);

% R_I_i=(rotz(rad2deg(Rd(3)))*roty(rad2deg(Rd(2)))*rotx(rad2deg(Rd(1))));
R_I_i=(roty(Rd(2))*rotx(Rd(1))*rotz(Rd(3)));
R_i_I = R_I_i';





%%
%���gantry�ĽǶ�
len=105;
xyz = [ r_xyz1;r_xyz2]; %��ʵ��ƽ������ϵI��XOY�����ж�Ӧ��ÿ�����������
stepAngle = 30;
step = 0:1:11;
ang = deg2rad(step' * stepAngle);
xw = radius*sin(ang);
yw = -radius*cos(ang);
zw = len*ones(size(xw));
xyzw = [[xw, yw, zw];[xw, yw, -zw]];  % ����������ϵW�еĸ���λ��
xyzi=((roty(Rd(2))*rotx(Rd(1)))*xyz')';
[gantry_angle,Zsi,Zdi]=optimization(xyzi,xyzw,Zsi,Zdi,Ysi,gantry_angle);
%%
Psi = [Xsi, Ysi, Zsi]';
PsI = R_i_I * Psi;
%���ƽ��̽������λ�ã�(Ydi,Zdi)���ο�����ϵΪƽ����������ϵ��Cho�����еģ�
%����Ford����II.C.4С�ڵ��������ƶ�(Y_d,Z_d)��ʾ����������ϵԭ����i����ϵ�е�����
%�ʣ�Xdi = Xwi��Ydi = Ywi�� Zdi = Zwi
Xwi = 0;
Zwi = Zdi;
Ywi = Ysi*Zdi/Zsi;
Pwi = [Xwi, Ywi, Zwi]';
%%
%ת��ƽ��̽������λ�úͷ������X��Դλ�õ���������ϵ��
%ת��ƽ��
%R_i_W  = rotaY(-90)*rotaX(90)*rotaZ(t); % tΪ����ƽ̨ת�ǣ���CBCT����ƽ̨��ת�Ƿ����෴
%t = rad2deg(angle); %����ƽ̨ת���Ƕ� ,�Ƕ���, TODO:ֱ�Ӽ���õ��ĽǶȲ�׼ȷ��
% t=rad2deg(gantry_angle-init_angle);
t=gantry_angle-init_angle;
 R_W_i = roty(-90)*rotx(90)*rotz(t);%��������ϵ������ƽ������ϵ��ת����
%R_W_i = rotz(90)*roty(-90);%��������ϵ������ƽ������ϵ��ת����
T_W_i = [R_W_i, Pwi; 0, 0, 0,1];%��������ϵ������ƽ������ϵ�任����
P_i_I = [0, 0, 0]';
T_i_I = [R_i_I, P_i_I; 0, 0, 0, 1]; %��ƽ����������ϵi��ƽ��ʵ������ϵI�ı任����
PsW = inv(T_W_i) * [Psi;1]; %X��Դ����������ϵW�е�λ��
% PsW=[roty(-t)*PsW(1:3);1];
T_I_W=inv(T_i_I*T_W_i);%��ʵƽ������ϵ����������ϵ��ת����

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
%��ƽ������ֱ��AB��CD�Ľ���
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
        disp('�޽��㣡');
        P = -1; 
else
        x=(C2*B1-C1*B2)/m;
        y=(C1*A2-C2*A1)/m;
        P = [x, y];
end

end

function [A,B,C] = GeneralEquation(first_x,first_y,second_x,second_y)
    %һ��ʽ Ax+By+C=0
    % from http://www.cnblogs.com/DHUtoBUAA/
    A=second_y-first_y;
    B=first_x-second_x;
    C=second_x*first_y-first_x*second_y;
end

function p = RotateAlongZ(point, angle)
%����������ת����ʵ��ƽ������ϵ��
%���룺   point, nx3,��n��1��x���꣬(n,2)y���꣬(n,3)z����
%         angleΪ��Z�����ת�Ƕ�,������
%�����   p Ϊnx3�ľ�����ת��ĵ�����
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
%����������ת����ʵ��ƽ������ϵ��
%���룺   pixel��i��1��ƽ�������u�����꣬pixel��i��2��λv�����꣬iΪ�к�
%        center(x,y)Ϊ��͸����������
%        scaleΪ���ش�С����λΪmm
%�����   p Ϊnx3�ľ���Ϊpixel��Ӧ������ʵ��ƽ������ϵ�е�����ֵ
[m, n] = size(pixel);
w=[];
for i=1:m
    % ��������ϵ��u��v����ƽ������ϵI��x,y,z��, u��x�����෴��v��y�ķ�����ͬ
    t = [-(pixel(i,1)-center(1))/scale, (pixel(i,2)-center(2))/scale];
    w = [w;t];
end
z = zeros(m,1);
p=[w,z];
end


function  ang = vectorAngleToX(v)
%�����������v��X��ļнǴ�С�뷽��
%   ���룺 v��                  1x3�������������
%
%   ����� ang��                �����������ƣ�v��X��ļнǣ���ʱ��Ϊ����ȡֵ��Χ[-pi/2,pi/2]��

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
    disp('��ע�⣬������������Ƕȣ�');
    eta = 0;
end

% n_z = [0,0,1];
% angle_d = acos(dot(c,n_z)/(norm(c)*norm(n_z))) % �жϽǶȷ�����ʱ��Ϊ����
% if ((angle_d >=0 ) && (angle_d < pi/2)) % ����
%     if ((0< angle) && (angle< pi/2))
%         eta = angle;
%     else
%         eta = angle-pi;
%     end
% else                                   %����
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

%�������
%ellipse equation: a*x^2 + b*x*y + c*y^2 +d*x + e*y + 1 = 0;
% Refer to��http://blog.csdn.net/ningyaliuhebei/article/details/46327681
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
K1= tan(theta1); %����б��
if (abs(K1)<1e-8) 
    K1=0;
end
%��Բ�����ĵ�
xc=(B*E-2*C*D)/(4*A*C-B^2);
yc=(B*D-2*A*E)/(4*A*C-B^2);
%��Բ�ĳ����̰��᳤��
la=sqrt(2*(A*xc^2+C*yc^2+B*xc*yc-1)/(A+C-sqrt((A-C)^2+B^2)));
lb=sqrt(2*(A*xc^2+C*yc^2+B*xc*yc-1)/(A+C+sqrt((A-C)^2+B^2)));
%�������ᶥ������
n1=[1,K1]; n1=n1/norm(n1);
center = [xc, yc];
va1=center+la*n1;  %���ᶥ��1, ����������
va2=center-la*n1;  %���ᶥ��2�������Ḻ��

if(K1 == 0)
    n2=[0,1];
else
    n2=[1,-1/K1];
end
n2=n2/norm(n2);
if n2(2) >= 0
    vb1=center+lb*n2;  %���ᶥ��1������������
    vb2=center-lb*n2;  %���ᶥ��2�������Ḻ��
else
    vb1=center-lb*n2;  %���ᶥ��1������������
    vb2=center+lb*n2;  %���ᶥ��2�������Ḻ��
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