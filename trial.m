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
syms u R w v v1 v2 w1 w2 real
use_simplified = false;

% define surface mapping function and its derivative
f = [R * cos(u) * cos(v), -R * cos(u) *sin(v), R*sin(u)];
fu = diff(f, u); fv = diff(f, v);
% define coordinate axes of normalized gauss frame
if(use_simplified)
    x = simplify(fu/norm(fu)); 
    y = simplify(fv/norm(fv)); 
    z = simplify(f/R);
else
    x = [-sin(u)*cos(v), sin(u)*sin(v), cos(u)].';
    y = [-sin(v), -cos(v), 0].';
    z = [cos(u) * cos(v), -cos(u)*sin(v), sin(u)];
end

% determine xh, yh 
xh = fu/R;
yh = fv/(R*cos(u));
zh = f/R;

% define normal vector to a point on the surface (here, a plane)
nh = f/R;

% define tensor forms for IAB and head
% Kh = zeros(2,2); Mh = eye(2,2); Th = zeros(1, 2);
% K_iab = [1/R, 0; 0, 1/R]; M_iab = [R 0; 0 R*cos(u)]; T_iab = [0, -tan(u)/R];

% curvatures from first principles (eqs 2.67 thesis)
norm_dfdu = R;
norm_dfdv = R*cos(u);

% metric forms
M_iab = diag([norm_dfdu, norm_dfdv]);
Mh = diag([1, 1]);

% curvature forms
K_iab = [xh.', yh.'].' * [(diff(nh, u)/norm_dfdu).', (diff(nh, v)/norm_dfdv).'];
K_iab = simplify(K_iab);
Kh = zeros(2,2);

%Torsion Equations
T_iab = yh * [(diff(xh, u)/norm_dfdu).', (diff(xh, v)/norm_dfdv).'];
T_iab = simplify(T_iab);
Th = zeros([1, 2]);

% define twist components
w = [w1, w2, 0].'; v = [v1, v2, 0].';
% w = [-w2, w1, 0].'; v = [v1, v2, 0].';

% determine rolling velocity of the head
fh = [0,0,0];
wt = simplify([xh.', yh.'].' * cross(nh, w).');
vt = simplify([xh.', yh.'].' * (cross(fh, w)+ v).');

% define relative rotational velocity projected to a surface normal
w_n = zh * w;

% solve it
syms p real;
%p = 0; % rolling on a plane; x axes coincide
Rpsi = [cos(p), -sin(p); -sin(p), -cos(p)];

alpha_o = simplify(alpha_head(Mh, Kh, K_iab, wt(1:2), vt(1:2), Rpsi))
alpha_f = simplify(alpha_iab(M_iab, Kh, K_iab, wt(1:2), vt(1:2), Rpsi))
psi_d = psi_dot(Th, Mh, T_iab, M_iab, alpha_o, alpha_f, w_n)