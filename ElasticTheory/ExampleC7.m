% MATLAB CODE: Example C-7
% Elasticity: Theory, Applications and Numerics 3e - Elsevier
% M.H. Sadd, University of Rhode Island
% Generates 2-D Warping Displacement Contours for 
% Elliptical Section Under Torsion - Figure 9-8
clc, clear all
% Input Geometry and Angle of Twist
a=1.0; b=0.5; alpha=1.0;
% Input Grid Space
[t,r]=meshgrid(0:pi/20:2*pi,0:0.05:1);
% Generate Contour Data
K=alpha*(b^2-a^2)/(a^2+b^2);
x=a*r.*cos(t);
y=b*r.*sin(t);
w=K*x.*y;
% Plotting Call with 20 Contours
contour(x,y,w,20,'k-');
axis equal