function alpha_iab = alpha_iab(Mi, Kh, K1, wt, vt)
% this function finds the local coordinates of the head
% Rpsi is the tangent ploane orientation
syms p real;
Rpsi = [cos(p), -sin(p); -sin(p), -cos(p)];
Miinv = pinv(Mi);
K1_tilde = Rpsi*K1*Rpsi;
rel_curve = pinv(Kh+K1_tilde);
alpha_iab = Miinv * Rpsi* rel_curve * (wt - Kh*vt);