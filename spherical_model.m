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

%% Multisphere Example I (Extension);
close all; clc; clear all
homedir = char(java.lang.System.getProperty('user.home'));
this_path = fullfile(homedir, 'Documents/screws');
savedir = 'Documents/Papers/SRS/Continuum/ICRA20/figures';
full_path=fullfile(homedir, savedir);

cd(this_path);
start = tic;
ri = 3.0/100; % meters
[Ri, Ro] = deal(2.7/100, 3/100); 

% material moduli for internal and external IAB skin
[C1, C2] = deal(1.1e4, 2.2e4); %eal(11000, 12000);
% IAB material density. Note pressure is already in psi
[rho, nu] = deal(38./386, .45);  % cue from dynamics of recatangular block example
[P_psi, model, ro, result] = bvp_free(C1, C2, Ri, Ro, rho, nu, 'extend', ri);
%P_psi = PressurePsi(P_Pa);
fprintf('Time to run: %f seconds', toc(start))
% Examine solution I
% Plot the displacements
figure(1)
subplot(131)
set(gca,'visible','off')
pdeplot3D(model,'ColorMapData',result.Displacement.ux);
title('x-displacement errors', 'FontSize', 30)
%colormap('jet')
colorbar('off')

subplot(132)
pdeplot3D(model,'ColorMapData',result.Displacement.uy);
xlabel('Time')
title('y-displacement errors', 'FontSize', 30)
%colormap('jet')
colorbar('off')
colorbar('south')

subplot(133)
pdeplot3D(model,'ColorMapData',result.Displacement.uz);
title('z-disp. errors', 'FontSize', 30)
%colormap('jet')
colorbar('off')


% subplot(131)
% pdegplot(model,'CellLabels','on','FaceAlpha',0.6,'FaceLabels','off')
% title('IAB Skins', 'FontSize', 30, 'FontWeight', 'bold');
% colorbar('off')

figure(2)
subplot(121)
pdeplot3D(model)
title('Mesh Model',  'FontSize', 30, 'FontWeight', 'bold');
colorbar('off')
colorbar('south')

subplot(122)
pdeplot3D(model,'ColorMapData',result.VonMisesStress);
title('Stress distribution', 'FontSize', 30, 'FontWeight', 'bold');
%colormap('jet')
colorbar('off')

%% Multisphere Example II (Compression);
close all; clc; 
homedir = char(java.lang.System.getProperty('user.home'));
this_path = fullfile(homedir, 'Documents/screws');

cd(this_path);

% Testing BVP for Contact free system
start = tic;
[Ri, Ro] = deal(2.5/100, 3/100);  % 50cm 70 cm
% material moduli for internal and external IAB skin
%[C1, C2] = deal(500000, 600000);
[C1, C2] = deal(1.1e4, 2.2e4);
% IAB material density
[rho, nu] = deal(4.8, .45);  % , .45 cue from dynamics of recatangular block example
ri = 2.25/100; % meters

[P_Pa, model, ro, result] = bvp_free(C1, C2, Ri, Ro, rho, nu, 'compress', ri);
P_psi = PressurePsi(P_Pa);
fprintf('Time to run: %f seconds', toc(start))
[ro_iso_chk, delta_volume, percent_volume] = isochoric_check(Ri, Ro, ri);
% Examine solution III
% Plot the displacements
figure(1)
subplot(131)
set(gca,'visible','off')
pdeplot3D(model,'ColorMapData',result.Displacement.ux);
title('x-disp. errors', 'FontSize', 30)
%colormap('jet')
colorbar('off')

subplot(132)
pdeplot3D(model,'ColorMapData',result.Displacement.uy);
xlabel('Time')
title('y-disp. errors', 'FontSize', 30)
%colormap('jet')
colorbar('off')
colorbar('south')

subplot(133)
pdeplot3D(model,'ColorMapData',result.Displacement.uz);
title('z-disp. errors', 'FontSize', 30)
%colormap('jet')
colorbar('off')


% subplot(131)
% pdegplot(model,'CellLabels','on','FaceAlpha',0.6,'FaceLabels','off')
% title('IAB Skins', 'FontSize', 30, 'FontWeight', 'bold');
% colorbar('off')

figure(2)
subplot(121)
pdeplot3D(model)
title('Mesh Model',  'FontSize', 30, 'FontWeight', 'bold');
colorbar('off')
colorbar('south')

subplot(122)
pdeplot3D(model,'ColorMapData',result.VonMisesStress);
title('Stress distribution', 'FontSize', 30, 'FontWeight', 'bold');
%colormap('jet')
colorbar('off')
