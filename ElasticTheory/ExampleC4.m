% MATLAB CODE: Example C-4
% Elasticity: Theory, Applications and Numerics 3e - Elsevier
% M.H. Sadd, University of Rhode Island
% Calculate and Plot Normalized Hoop Stress on Circular Hole
% In Infinite Plane Under Uniform Tension at Infinity: Example 8.7
% Non-Dimensional Plot Generates Figure 8-13
clc;clear all;
% Input (r/a)- Variable and Generate Angular Coordinate Space
r=1;
t=[0:0.01:2*pi];
% Calculation Loop
st=0.5*(1+(1/r)^2)-0.5*(1+3*(1/r)^4)*cos(2*t);
% Plotting Call
polar(t,st)
title('Non-Dimensional Hoop-stress Around Hole')