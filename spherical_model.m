%% dyadics in spherical coordinates
clear all; clc
syms ri Ri R T Ph ph t p C1 C2  real
r = (R^3 + ri^3 - Ri^3)^(1/3); t = T; ph = Ph;
F = [R^2/r^2   0     0;...
           0         r/R     0; ...
           0            0           r/R];  
B = F*F.';
C = F.'*F;

StressMat = -p * eye(3) + C1*C + C2 * pinv(C);
sigma_rr = StressMat(1,1);
sigma_pp = StressMat(2,2);
sigma_tt = StressMat(3,3);

% confirm \sigma_{rr} + sigma_pp - 2 sigma_tt
integ_int = sigma_tt + sigma_pp - 2 * sigma_rr;
simplify(integ_int)