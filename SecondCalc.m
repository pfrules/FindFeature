function Parameter=SecondCalc(FirstOutput,pixel_scale)

for i=1:size(FirstOutput,1)
PS(:,i)=FirstOutput(i).PsW; 
TD(:,i)=FirstOutput(i).T_I_W(:,4);
end
IsocentricPoint=mean(PS,2);
[normalz,p]=PlaneDirectFit(PS(1:3,:)');
SAD=mean(sqrt(sum((PS(1:3,:)-repmat(IsocentricPoint(1:3),1,size(PS,2))).^2,1)));
nx=[p;1]-IsocentricPoint;
normalx=(nx(1:3)/norm(nx))';
normaly=-cross(normalx,normalz);
R_W0_W=[normalx;normaly;normalz];
[U,D,V]=svd(R_W0_W);
R_W0_W=U*V';
T_W_W0=[[R_W0_W';0 0 0],IsocentricPoint];%等中心点世界坐标系转换到体膜坐标系变换矩阵
% Psw=T_W_W0*PS;

% ooo2=mean(TD,2);
% hold on
% plot3(PS(1,:),PS(2,:),PS(3,:),'b.','MarkerSize',0.5) 
% plot3(IsocentricPoint(1,1),IsocentricPoint(2,1),IsocentricPoint(3,1),'b.','MarkerSize',10)
% plot3(TD(1,:),TD(2,:),TD(3,:),'r.','MarkerSize',0.5) 
% plot3(ooo2(1,1),ooo2(2,1),ooo2(3,1),'r.','MarkerSize',10) 
%%
Parameter=[];
for i=1:size(FirstOutput,1)
    T_W_i=FirstOutput(i).T_W_i;
    IsocentricPoint_i=T_W_i*IsocentricPoint;
    Psi=FirstOutput(i).Psi;
    R_i_I=FirstOutput(i).R_i_I;                         %原虚拟平板到原现实平板旋转矩阵
    IsocentricPoint_I=getpoint(R_i_I',IsocentricPoint_i(1:3),Psi);
    Pixel=FirstOutput(i).Pixel;
    center=[Pixel(1)-IsocentricPoint_I(1)*pixel_scale;Pixel(2)+IsocentricPoint_I(2)*pixel_scale];
    T_I0_W0=FirstOutput(i).T_I_W;                       %原现实平板到原世界坐标系变换矩阵
    T_I_I0=[eye(3),[IsocentricPoint_I';0];[0 0 0 1]];   %现实平板到原现实平板变换矩阵
    T_I_W=inv(T_W_W0)*T_I0_W0*T_I_I0;                   %现实平板到世界坐标系变换矩阵

    PsI=inv(T_I_I0)*[FirstOutput(i).PsI;1];
    PsW=T_I_W*PsI;
    t=acos(PsW(1)*SAD/(SAD*norm(PsW(1:2))))*180/pi;
    if abs(FirstOutput(i).t)>180
        t=360-t;
    end
%     RiI=(((rotz(90)*roty(-90)*rotz(-t))*T_I_W(1:3,1:3)))';
    R_W_i=roty(-90)*rotx(90)*(rotz(sign(FirstOutput(i).t)*t));
     RiI=(R_W_i*T_I_W(1:3,1:3))';
%      RiI= ((rotz(90)*roty(-90)*rotz(FirstOutput(i).t))*T_I_W(1:3,1:3));
    theta=asin(-RiI(3,2)/sqrt(RiI(3,1)^2+RiI(3,3)^2))*180/pi;
    phi=atan(RiI(3,1)/RiI(3,3))*180/pi;
    eta=atan(RiI(1,2)/RiI(2,2))*180/pi;
    %PsW=inv(T_W_W0)*FirstOutput(i).PsW;
    parameter.num=FirstOutput(i).num;
    parameter.Pixel=center;                             %等中心点对应的投影点
%   parameter.Rd=FirstOutput(i).Rd;                     %平板相对体膜世界坐标系的偏转角度
    parameter.T_I_W=T_I_W;                              %现实平板到世界坐标系变换矩阵
    parameter.PsI=PsI;                                  %现实平板下源的坐标系
    parameter.PsW=PsW;                                  %世界坐标系下源的坐标
    parameter.DW=T_I_W*[0;0;0;1];
    parameter.SAD=SAD;                                  %源轴距
    parameter.SID=PsI(3);                               %源像距
    parameter.theta=theta;
    parameter.phi=phi;
    parameter.eta=eta;
    parameter.t=t;
    Parameter = [Parameter; parameter];
	clc;
    disp(i);
end 

end
     

function p=getpoint(R,x1,x2)
parameter=R(:,3);
P12=x2-x1;
n=-(parameter(1)*(x2(1)-x1(1))+parameter(2)*(x2(2)-x1(2))+parameter(3)*(x2(3)-x1(3)))/(parameter(1)*x1(1)+parameter(2)*x1(2)+parameter(3)*x1(3));
p3=x1+P12/n;
p3=R'*p3;
p=p3(1:2)';
end

function [normal,p]=PlaneDirectFit(planeData)
%空间平面拟合
xyz0=mean(planeData,1);
centeredPlane=bsxfun(@minus,planeData,xyz0);
[U,S,V]=svd(centeredPlane);
a=V(1,3);
b=V(2,3);
c=V(3,3);
d=-dot([a b c],xyz0);
normal=[a,b,c];
if c<0
    normal=-normal;
end
t=[a,b,c,d]*[planeData(1,:)';1];
p=planeData(1,:)'-[a*t;b*t;c*t];
end

