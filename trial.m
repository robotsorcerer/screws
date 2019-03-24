%% Page 88 of Murray and Sastry
% Computes example exponential coordinates for twists
clc, clear all
[w1 w2 w3] = deal([0 0 1].');

syms l1 l2 real;
% points on axes wi's
q1 = [0 0 0]';
q2 = [0 l1, 0]';
q3 = [0 l1+l2 0]';
% twists
z1 = twist(w1, q1);
z2 = twist(w2, q2);
z3 = twist(w3, q3);

[v1 v2 v3] = deal(z1(1:3), z2(1:3), z3(1:3));

% exponential maps for twists
syms theta1 theta2 theta3 theta4 v4 real
ez1 = emap_twist(w1, v1, theta1);
ez2 = emap_twist(w2, v2, theta2);
ez3 = emap_twist(w3, v3, theta3);

w4 = [0 0 0]'; 
v4 = [0 0 1]';
ez4 = emap_twist(w4, v4, theta4);

R = ez1(1:3,1:3) * ez2(1:3, 1:3) * ez3(1:3,1:3);
R = simplify(R);

p = ez1(1:3,4) + ez2(1:3, 4) + ez3(1:3,4);
p = simplify(p);

% transformation vector
g123 = ez1*ez2*ez3*ez4;
g123 = simplify(g123)

%% Planar twist
clc
syms qx qy r l1 l2 h t11 t12 t21 t22 t real
w = [0 0 1]';
q = [qx qy 0]';

z00 = twist(w, subs(q, [qx, qy], [0, h+l2]));
z11 = [0 r 0 0 0 1]'; z12 = [l1 r 0 0 0 1]'; 
z21 = [h -r 0 0 0 1]'; z22 = [h+l2; -r; 0; 0; 0; 1];

gst1 = emap_twist(z11, t11) * emap_twist(z12, t12); % ...
gst2 = emap_twist(z21, t21) * emap_twist(z22, t22);

gst0 = emap_twist(z00, t)
gst1 = simplify(gst1)
gst2 = simplify(gst2)

%% Analysis of deformation: Audrey Sedal
clc; clear all
syms r ri R Ri  ell L t T z Z
lambdaz = ell/L;

r = sqrt( (R^2 - Ri^2)/(ell/L) + ri^2 );
z = lambdaz*Z;
t = T;

F = [diff(r, R), diff(r, T) diff(r, Z); ...
     diff(t, R), diff(t, T) diff(t, Z);
     diff(z, R), diff(z, T) diff(z, Z)];
    

lambda1 = R*L/(r*ell);
lambda2 = r/R;
lambda3 = ell/L;

I2 = simplify(1/lambda1^2 + 1/lambda2^2 + lambda1^2*lambda2^2);
