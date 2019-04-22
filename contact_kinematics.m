%% General computaions from Shastry and Murray's 1990 paper 
% on multifingered robot hands
clc; clear all
syms r p 
G = [eye(3) zeros(3,3); skewsem([r*cos(p), r*sin(p), 0]), eye(3)];
F = zeros(6, 2);
F(1,1) = 1;
F(2, 2) = 1;

Fo = G * F