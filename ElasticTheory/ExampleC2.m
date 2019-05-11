% MATLAB CODE: Example C-2
% Elasticity: Theory, Applications and Numerics 3e - Elsevier
% M.H. Sadd, University of Rhode Island
% Program to Numerically Calculate Invariants, Principal Values
% and Directions of a Matrix 
% Program Uses Matrix from Example 1-3
clc;clear all;
% Input Matrix
A=[2,0,0;0,3,4;0,4,-3]
% Calculate Invariants
invariants=[trace(A),(trace(A)^2-trace(A*A))/2,det(A)]
[V,L]=eig(A);
% Principal Values are the Diagonal Elements of the L Matix
principal_values=[L(1,1),L(2,2),L(3,3)]
% Principal Directions are the Columns of the V Matrix
principal_directions=[V(:,1),V(:,2),V(:,3)]