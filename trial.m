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

%% Spherical example
clc, clear all
syms R u v real
f = [R * cos(u) * cos(v), -R * cos(u) *sin(v), R*sin(u)];
fu = diff(f, u); fv = diff(f, v);
x = simplify(fu/norm(fu)); y = simplify(fv/norm(fv)); z = simplify(f/R);
% Curvature form
% K = [x(u), y(u)]^T * [dz(u)/du)/norm($\partial f/\partial u$),
%                       dz(u)/du)/norm($\partial f/\partial v$)]$
% K = curvature_form(f);
x = [-sin(u)*cos(v), sin(u)*sin(v), cos(u)].';
y = [-sin(v), -cos(v), 0].';
z = [cos(u) * cos(v), -cos(u)*sin(v), sin(u)];
K = simplify([x, y].' * [(diff(z, u)/norm(fu)).', (diff(z, u)/norm(fv)).']);
T = simplify(y.' * [diff(x, u)/norm(fu), diff(x, v)/norm(fv)]);

syms p real;
Rpsi = [cos(p), -sin(p); -sin(p), -cos(p)];
%% Example Murray book, 1990
clc; clear all
syms u R v1 v2 w1 w2 real
Ko = zeros(2,2); Mo = eye(2,2); To = zeros(1, 2);
Kf = [1/R, 0; 0, 1/R]; Mf = [R 0; 0 R*cos(u)]; Tf = [0, -tan(u)/R];

wt = [w1, w2].'; vt = [v1, v2].';
alpha_o = simplify(alpha_head(Mo, Ko, Kf, wt, vt))
alpha_f = simplify(alpha_iab(Mf, Ko, Kf, wt, vt))
psi_d = psi_dot(To, Mo, Tf, Mf, alpha_o, alpha_f)

% syms p real;
% Rpsi = [cos(p), -sin(p); -sin(p), -cos(p)];
% Kf_tilde = Rpsi*Kf*Rpsi;
% rel_curve = pinv(Ko+Kf_tilde);
% alpha_h = pinv(Mo) * rel_curve * (wt - Kf_tilde*vt);

%% 
clc; clear all
syms r p real
B = zeros(6, 2);
B(1,1) = 1; B(2, 2) = 1;
G = [eye(3), zeros(3,3); skewsem([r*cos(p), r*sin(p), 0]), eye(3)] * B

%% IAB dynamics
clc; clear all
syms r l t p
rr = [r*cos(t)*sin(p); r*sin(t)*sin(p); r*cos(p)];
% rrd = simplify(diff(rr, r) + diff(rr, t)+ diff(rr, p));
syms rd td pd
rrd =  [rd*cos(t)*sin(p) + r*pd*cos(p)*cos(t) - r*td*sin(p)*sin(t);
 rd*sin(p)*sin(t) + r*pd*cos(p)*sin(t) + r*td*cos(t)*sin(p);
                                 rd*cos(p) - r*pd*sin(p)]
                             
 %% 
syms r th p rho
syms r(t)
diff(rho * diff(r))
pretty(ans)
