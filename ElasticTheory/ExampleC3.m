% MATLAB CODE: Example C-3
% Elasticity: Theory, Applications and Numerics 3e - Elsevier
% M.H. Sadd, University of Rhode Island
% Calculate and Plot Stresses in Thick-Walled Cylinder Problem Example 8.6
% Internal Pressure Loading Case with r1/r2=0.5; Generate Figure 8-9
clc;clear all;
% Generate Non-Dimensional Radial Coordinate Space: r/r2
r=[0.5:0.01:1];
% Calculation Loop for Stresses
for i=1:length(r)
    sr(i)=(1/3)*(1-(1/r(i))^2);
    st(i)=(1/3)*(1+(1/r(i))^2);
end
% Plotting Call
plot(r,sr,r,st)
xlabel('Dimensionless Distance, r/r_2')
ylabel('Dimensionless Stress')
