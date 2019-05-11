% MATLAB CODE: Example C-6
% Elasticity: Theory, Applications and Numerics 3e - Elsevier
% M.H. Sadd, University of Rhode Island
% Calculate and Plot Tau Stress Contours Under Flat Punch Loading (P=1)
% Numerically Evaluate Integrals in Solution (8.5.9)Using quadv(.)
% Singularity at s=a is Handled using Integral Range Limiter e
function flat_punch
clc;clear all;clf
a=1; %input loading width
e=0.0001; %input range limiter
[x,y]=meshgrid(-2*a:0.1:2*a,0.05:0.1:4*a+0.1);
for i=1:length(x)
    for j=1:length(y)
        sx(i,j)=quadv(@(s)Irx(x(i,j),y(i,j),s,a),-(a-e),a-e);
        sy(i,j)=quadv(@(s)Iry(x(i,j),y(i,j),s,a),-(a-e),a-e);
        txy(i,j)=quadv(@(s)Irxy(x(i,j),y(i,j),s,a),-(a-e),a-e);
        smax(i,j)=sqrt((((sx(i,j)-sy(i,j))/2)^2)+(txy(i,j)^2));
    end
end
% Draw half-space boundary line
plot([-2*a,2*a],[0,0],'k','linewidth',3)
hold on
contour(x,-y,smax,40,'k','linewidth',1.4)
axis off;hold off;
title('\tau_m_a_x Contours: Frictionless Rigid Punch Loading on a Half-Space')
function I=Irx(x,y,s,a)
I=-(2*y/pi)*(1/(pi*sqrt(a^2-s^2)))*(x-s)^2/(((x-s)^2+y^2)^2);
function I=Iry(x,y,s,a)
I=-(2*y^3/pi)*(1/(pi*sqrt(a^2-s^2)))/(((x-s)^2+y^2)^2);
function I=Irxy(x,y,s,a)
I=-(2*y^2/pi)*(1/(pi*sqrt(a^2-s^2)))*(x-s)/(((x-s)^2+y^2)^2);
