function alpha_h = alpha_head(Mh, Kh, K1, wt, vt, Rpsi)
% this function finds the local coordinates of the head
% Rpsi is the tangent ploane orientation
K1_tilde = Rpsi*K1*Rpsi;
rel_curve = pinv(Kh+K1_tilde);
alpha_h = pinv(Mh) * rel_curve * (wt - K1_tilde*vt);
%alpha_h = pinv(Mh) * rel_curve * (wt);