function alpha_iab = alpha_iab(Mi, Kh, K1, wt, vt, Rpsi)
% this function finds the local coordinates of the head
% Rpsi is the tangent ploane orientation
% syms p real;
% Rpsi = [cos(p), -sin(p); -sin(p), -cos(p)];
K1_tilde = Rpsi*K1*Rpsi;
rel_curve = pinv(Kh+K1_tilde);
alpha_iab = pinv(Mi) * Rpsi* rel_curve * (wt - Kh*vt);
%alpha_iab = pinv(Mi) * Rpsi* rel_curve * (wt);