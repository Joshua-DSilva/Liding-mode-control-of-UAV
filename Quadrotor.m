function dx = QRBS(t,x)
m = 486e-3;d = 25e-2;
ob = x(13)-x(14)+x(15)-x(16);
beta0 = 189.63; beta1 = 6.0612; beta2 = 0.0122; b = 280.19;
g = 9.8; c = 1;
Ix = 3.8278e-3; Iy = 3.8288e-3; Iz = 7.6566e-3;
Kfax = 5.567e-4;Kfay = 5.567e-4;Kfaz = 6.354e-4;
Kftx = 5.567e-4;Kfty = 5.567e-4;Kftz = 6.354e-4;
Jr = 2.8385e-5;
a1 = (Iy-Iz)/Ix;
a2 = -Kfax/Ix;
a3 = -Jr/Ix;
a4 = (Iz-Ix)/Iy;
a5 = -Kfay/Iy;
a6 = Jr/Iy;
a7 = (Ix-Iy)/Iz;
a8 = -Kfaz/Iz;
a9 = -Kftx/m;
a10 = -Kfty/m;
a11 = -Kftz/m;
b1 = d/Ix;
b2 = d/Iy;
b3 = 1/Iz;
xd = [0 0 0 0 0 0 sin(t) cos(t) cos(t) -sin(t) 0.1*t 0.1]';
% [0 0 0 0 0 0 sin(t) cos(t) cos(t) -sin(t) 0.1*t 0.1]';
z = zeros(12,1); xdd = [0 0 0 0 0 0 cos(t) -sin(t) -sin(t) -cos(t) 0.1 0]';
% z = zeros(12,1); xdd = [0 0 0 0 0 0 cos(t) -sin(t) -sin(t) -cos(t) 0.1 0]';
for w = 1:2:11
    z(w) = xd(w) - x(w);
end
for q = 2:2:12
    z(q) = x(q) - xdd(q-1) - z(q-1);
end
Cp = 1;
Ct = 1;
V1 = z(1)^2/2;
V2 = (V1 + z(2)^2)/2;
V3 = z(3)^2/2;
V4 = (V3 + z(4)^2)/2;
U2 = 1/b1*(-sign(z(2))-z(2)-a1*x(4)*x(6)-a2*x(2)^2-a3*ob*x(4)+xdd(2)+(xd(2)-x(2)));
U3 = 1/b2*(-sign(z(4))-z(4)-a4*x(2)*x(6)-a5*x(4)^2-a6*ob*x(2)+xdd(4)+(xd(4)-x(4)));
U4 = 1/b3*(-sign(z(6))-z(6)-a7*x(2)*x(4)-a8*x(6)^2+xdd(6)+(xd(6)-x(6)));
U1 = m/(Cp*Ct)*(-sign(z(12))-z(12)-a11*x(12)+xdd(12)+(xd(12)-x(12))+g);
if U1==0
    Ux = 0;
    Uy = 0;
else
    Ux = m/U1*(-sign(z(8))-z(8)-a9*x(8)+xdd(8)+(xd(8)-x(8)));
    Uy = m/U1*(-sign(z(10))-z(10)-a10*x(10)+xdd(10)+(xd(10)-x(10)));
end
x1dot = x(2);
x2dot = a1*x(4)*x(6)+a2*x(2)^2+a3*ob*x(4)+b1*U2;
x3dot = x(4);
x4dot = a4*x(2)*x(6)+a5*x(4)^2+a6*ob*x(2)+b2*U3;
x5dot = x(6);
x6dot = a7*x(2)*x(4)+a8*x(6)^2+b3*U4;
x7dot = x(8);
x8dot = a9*x(8) + Ux*U1/m;
x9dot = x(10);
x10dot = a10*x(10) + Uy*U1/m;
x11dot = x(12);
x12dot = a11*x(12)+ Ct*Cp*U1/m-g;
x13dot = b*V1-beta0-beta1*x(13)-beta2*x(13)^2;
x14dot = b*V2-beta0-beta1*x(14)-beta2*x(14)^2;
x15dot = b*V3-beta0-beta1*x(15)-beta2*x(15)^2;
x16dot = b*V4-beta0-beta1*x(16)-beta2*x(16)^2;
dx = [x1dot;x2dot;x3dot;x4dot;x5dot;x6dot;x7dot;x8dot;x9dot;x10dot;x11dot;x12dot;x13dot;x14dot;x15dot;x16dot];



