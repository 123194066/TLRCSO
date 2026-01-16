clc;clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Jx与Jy推导的公式，Jx即功率对x的偏导数
syms x y z x_t y_t z_t Pt  miu_T rho_1 
lamba_U=0.3279;G_T =1.58;G_R = 3.98;
% gamma_1 = -sqrt(0.1);tao = 0.5;
X=x_t-x;
Y=y_t-y;
Z=z-z_t;
d0=sqrt(X.^2+Y.^2+Z.^2);
% d1 = d0./ cos(atan(2 * z / d0));
zr=0;  %%%  x,y,z为与地面的反射点
tt=(zr-z_t)./(-z-z_t);
coordinate(1)=tt.*(x-x_t)+x_t;
coordinate(2)=tt.*(y-y_t)+y_t;
d1=sqrt((coordinate(1)-x).^2+(coordinate(2)-y).^2+(0-z).^2); %阅读器和地面反射点之间的距离
d2=sqrt((x_t-coordinate(1)).^2+(y_t-coordinate(2)).^2+(z_t-0).^2);%标签和地面反射点之间的距离
 Ld01=1./(d0).^2.*((1+0.3.*d0./(d1+d2).*cos(2.*pi.*(d1+d2-d0)./lamba_U)).^2+(0.3 .*d0./(d1+d2).*sin(2.*pi.*(d1+d2-d0)./lamba_U)).^2);
 P11=0.25.* Pt.*G_T.^2.*G_R.^2 .*(abs(Ld01)).^2;
Jx=diff(P11,x_t)
Jy=diff(P11,y_t)