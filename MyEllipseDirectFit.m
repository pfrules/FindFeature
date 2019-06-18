function [Re, center, vertex]= MyEllipseDirectFit(XY)
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

Re=A;

%长轴倾角
%ellipse equation: a*x^2 + b*x*y + c*y^2 +d*x + e*y + 1 = 0;
% Refer to：http://blog.csdn.net/ningyaliuhebei/article/details/46327681
B=Re(2)/Re(6); A=Re(1)/Re(6);C=Re(3)/Re(6); D=Re(4)/Re(6);E=Re(5)/Re(6);
theta1=0.5*atan(B/(A-C)); 
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