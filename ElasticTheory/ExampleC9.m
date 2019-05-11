% MATLAB CODE: Example C-9
% Elasticity: Theory, Applications and Numerics 3e - Elsevier
% M.H. Sadd, University of Rhode Island
% Calculate Beta Parameters for Orthotropic E-Glass Material
clc; clear all;
% Input Number of Materials, Names and Stiffness Moduli
N=1;
name(1,:)='E-Glass/Epoxy ';
e1(1)=38.6; e2(1)=8.3; nu12(1)=0.26; mu12(1)=4.2;
% Calulate Compliances
for i=1:N
s11(i)=1/e1(i); s22(i)=1/e2(i); s12(i)=-nu12(i)/e1(i); s66(i)=1/mu12(i);
% Calculate Beta Values
b1(i)=sqrt(-(1/(2*s11(i)))*(-(2*s12(i)+s66(i))+sqrt((2*s12(i)+s66(i))^2-4*s11(i)*s22(i))));
b2(i)=sqrt(-(1/(2*s11(i)))*(-(2*s12(i)+s66(i))-sqrt((2*s12(i)+s66(i))^2-4*s11(i)*s22(i))));
% Print Results to Screen
fprintf(1,'\n ')
disp(name(i,:))
fprintf(1,'    beta(1)=%5.3f',b1(i))
fprintf(1,'    beta(2)=%5.3f',b2(i))
end