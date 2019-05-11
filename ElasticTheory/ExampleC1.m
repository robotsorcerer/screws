% MATLAB CODE: Example C-1
% Elasticity: Theory, Applications and Numerics 3e - Elsevier
% M.H. Sadd, University of Rhode Island
% Program to Calculate Components of Matrix (Second Order Tensor) Under Rotational Transformation
% Q = Rotation Matrix, A = Original Matrix to be Transformed, AP = Transformed Matrix
clc;clear all;
% Input Matrices in MATLAB Format
Q=[0,0,-1;0,1,0;1,0,0];
A=[1,1,1;0,4,2;0,1,1];
% Apply Transformation Law
AP=Q*A*Q';
% Display Segement
disp('Original Matrix')
disp(A)
disp('Rotation Matrix')
disp(Q)
disp('Transformed Matrix')
disp(AP)
