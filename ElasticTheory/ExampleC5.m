% MATLAB CODE: Example C-5
% Elasticity: Theory, Applications and Numerics 3e - Elsevier
% M.H. Sadd, University of Rhode Island
% Displacement Vector Distribution Plot for Flamant Problem - Figure 8-19
% Normal Loading Case (X=0), Generates Figure 8-22
clc;clear all;
% Input Parameters
Y=1; E=1; nu=0.3;
% Input Radial Coordinates: logspace(M,N,*) Generates Region 10^M < r < 10^N
r=logspace(-3,-0.5,40);
% Input Theta Coordinates
t=[0:pi/20:pi];
% Create Mesh
[t,r]=meshgrid(t,r);
ur=Y/(pi*E)*((1-nu)*(t-pi/2).*cos(t)-2*log(r).*sin(t));
ut=Y/(pi*E)*(-(1-nu)*(t-pi/2).*sin(t)-2*log(r).*cos(t)-(1+nu)*cos(t));
% Calculate Cartesian Displacement Components - Flip y-Component
ux=cos(t).*ur-sin(t).*ut;
uy=-(sin(t).*ur+cos(t).*ut);
% Covert to Cartesian Coordinate Mesh
[x,y]=pol2cart(-t,r);
% Plotting Call for Vector Distribution
quiver(x,y,ux,uy)
axis equal;