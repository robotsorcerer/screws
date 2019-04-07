function dot_psi = psi_dot(Th, Mh, T_iab, M_iab, alpha_head_dot, alpha_iab_dot)
% this function finds the angles between t6hye two local coordinates
% Rpsi is the tangent ploane orientation
dot_psi = Th*Mh*alpha_head_dot + T_iab*M_iab*alpha_iab_dot