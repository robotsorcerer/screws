% MATLAB CODE: Example C-8
% Elasticity: Theory, Applications & Numerics 3e - Elsevier
% M.H. Sadd, University of Rhode Island
% Three Dimensional Plot of Warping Displacement Surface
% for Elliptical Section Under Torsion
clc;clear all;
a=1;b=0.5;
[t,r]=meshgrid(0:pi/20:2*pi,0:0.05:1);
x=a*r.*cos(t);
y=b*r.*sin(t);
w=-x.*y;
surfc(x,y,w);
colormap gray
h=findobj(gcf,'type','patch');
set(h,'LineWidth',1.0, 'EdgeColor','k');
axis([-1 1 -1 1])
axis square
view(20,10)