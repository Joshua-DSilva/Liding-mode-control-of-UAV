clc;clear all;close all;
[t,x] = ode45(@Quadrotor,[0,15],[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
plot3(x(:,7),x(:,9),x(:,11),':b',LineWidth=1.5)
hold on
plot3(sin(t),cos(t),0.1*t,':r',LineWidth=1.5)
%plot3(cos(t),sin(t),1+0*t,':r',LineWidth=1.5)
% legend([plt1(1),plt2(1)])