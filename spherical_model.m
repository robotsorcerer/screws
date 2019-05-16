%% dyadics in spherical coordinates
clear all; clc
syms ri Ri R T Ph ph t p C1 C2  real
r = (Ro^3 + ri^3 - Ri^3)^(1/3); t = T; ph = Ph;
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

%% Deformation analysis
% some useful examples
% https://www.mathworks.com/help/pde/examples/dynamic-analysis-of-a-clamped-beam.html
clc; clear all; close all
%profile on -history

homedir = '/home/lex/';
otherpath = 'Documents/superchicko/sofa/IAB8/compoz';
sphere_path=fullfile(homedir, otherpath, 'sphere.STL');
model = createpde(1);
importGeometry(model, sphere_path);

% Plots
close all
figure(1)
pdegplot(model,'FaceLabels','on')
view(30,30);
title('Solidworks Model')

applyBoundaryCondition(model,'dirichlet','face',1,'u',0);

% material moduli for internal and external IAB skin
[C1, C2] = deal(11000, 22000);
[rho, nu] = deal(.5/386., .45); 
% specify material moduli, poison ratio for the deformation;
fun = @(R, r) 2.* C1*((1./r)-((R.^6)./(r.^7))) + ...
           2.* C2*((r.^5)./(R.^6) - (R.^10)./(r.^11));

P = integral(@(R)fun(R, ro), Ri, Ro);
    
% define the boundary conditions
structuralBC(model,'Face',1,'Constraint','symmetric'); % symmetric isochoric

%Specify the gravity load on the sphere.
structuralBodyLoad(model,'GravitationalAcceleration',[0;0;-9.8]);

%% Multisphere Example I (Extension);
close all; clc; clear all
homedir = '/home/lex/';
savedir = 'Documents/Papers/PhDThesis/figures/deformation';

% Testing BVP for Contact free system
cd('/home/lex/Documents/MATLAB/screws');
start = tic;
[Ri, Ro] = deal(.1, .15);  % 50cm 70 cm

% material moduli for internal and external IAB skin
[C1, C2] = deal(500000, 600000); %eal(11000, 12000);
% IAB material density
[rho, nu] = deal(.38/386., .487);  % cue from dynamics of recatangular block example
ri = .13; % meters
[P_Pa, model, ro, result] = bvp_free(C1, C2, Ri, Ro, rho, nu, 'extend', ri);
P_psi = PressurePsi(P_Pa);
fprintf('Time to run: %f seconds', toc(start))
% Examine solution I
%homedir = '/Users/olalekanogunmolu/';
homedir = '/home/lex/';
full_path=fullfile(homedir, savedir);
% Plot the displacements
figure(1)
subplot(131)
pdeplot3D(model,'ColorMapData',result.Displacement.ux*100);
title('x-displacement')
% : R_i = 50cm, R_o = 70cm, r_i = R_i + 20 (cm)
colormap('jet')
% colorbar('off')

subplot(132)
pdeplot3D(model,'ColorMapData',result.Displacement.uy*100);
xlabel('Time')
title('y-displacement')
colormap('jet')
% colorbar('off')

subplot(133)
pdeplot3D(model,'ColorMapData',result.Displacement.uz*100);
title('z-displacement')
% : R_i = 50cm, R_o = 70cm, r_i = R_i + 20 (cm)
colormap('jet')

figure(2)
subplot(131)
pdegplot(model,'CellLabels','on','FaceAlpha',0.6,'FaceLabels','off')
title('IAB Skins');

subplot(132)
pdeplot3D(model)
title('Mesh Model')

subplot(133)
pdeplot3D(model,'ColorMapData',result.VonMisesStress);
title('Stress distribution')
colormap('jet')

%% Multisphere example II (Extension)
close all; clc; 

start = tic;
[Ri, Ro] = deal(.075, .1);  % 50cm 70 cm

% material moduli for internal and external IAB skin
[C1, C2] =  deal(500000, 600000);
% IAB material density
[rho, nu] = deal(1e-4, .45);  %deal(.38/386., .4995); %% cue from dynamics of recatangular block example
ri = .12; % meters
close all
[P, model, ro, result] = bvp_free(C1, C2, Ri, Ro, rho, nu, 'extend', ri)
fprintf('Time to run: %f seconds', toc(start));

homedir = '/home/lex/';
% homedir = '/Users/olalekanogunmolu/';
full_path=fullfile(homedir, savedir);
% Plot the displacements
figure(1)
subplot(131)
pdeplot3D(model,'ColorMapData',result.Displacement.ux*100);
title('x-displacement')
% : R_i = 50cm, R_o = 70cm, r_i = R_i + 20 (cm)
colormap('jet')
colorbar('off')

subplot(132)
pdeplot3D(model,'ColorMapData',result.Displacement.uy*100);
xlabel('Time')
title('y-displacement')
colormap('jet')
colorbar('off')

subplot(133)
pdeplot3D(model,'ColorMapData',result.Displacement.uz*100);
title('z-displacement')
% : R_i = 50cm, R_o = 70cm, r_i = R_i + 20 (cm)
colormap('jet')

figure(2)
subplot(131)
pdegplot(model,'CellLabels','on','FaceAlpha',0.6,'FaceLabels','off')
title('IAB Skins');

subplot(132)
pdeplot3D(model)
title('Mesh Model')

subplot(133)
pdeplot3D(model,'ColorMapData',result.VonMisesStress);
title('Stress distribution')
colormap('jet')

%% Multisphere Example III (Compression);
close all; clc; clear all
homedir = '/home/lex/';
%homedir = '/Users/olalekanogunmolu/';
savedir = 'Documents/Papers/PhDThesis/figures/deformation';

% Testing BVP for Contact free system
cd('/home/lex/Documents/MATLAB/screws');
start = tic;
[Ri, Ro] = deal(.12, .15);  % 50cm 70 cm

% material moduli for internal and external IAB skin
[C1, C2] = deal(500000, 1200000);
% IAB material density
[rho, nu] = deal(.85/386000, .45);  % , .45 cue from dynamics of recatangular block example
ri = .1; % meters
close all
[P_Pa, model, ro, result] = bvp_free(C1, C2, Ri, Ro, rho, nu, 'compress', ri);
P_psi = PressurePsi(P_Pa);
fprintf('Time to run: %f seconds', toc(start))
% Examine solution III
homedir = '/home/lex/';
full_path=fullfile(homedir, savedir);
% Plot the displacements
figure(1)
subplot(131)
pdeplot3D(model,'ColorMapData',result.Displacement.ux*1000);
title('x-displacement')
% : R_i = 50cm, R_o = 70cm, r_i = R_i + 20 (cm)
colormap('jet')
% colorbar('off')

subplot(132)
pdeplot3D(model,'ColorMapData',result.Displacement.uy*1000);
xlabel('Time')
title('y-displacement')
colormap('jet')
% colorbar('off')

subplot(133)
pdeplot3D(model,'ColorMapData',result.Displacement.uz*1000);
title('z-displacement')
% : R_i = 50cm, R_o = 70cm, r_i = R_i + 20 (cm)
colormap('jet')

figure(2)
subplot(131)
pdegplot(model,'CellLabels','on','FaceAlpha',0.6,'FaceLabels','off')
title('IAB Skins');

subplot(132)
pdeplot3D(model)
title('Mesh Model')

subplot(133)
pdeplot3D(model,'ColorMapData',result.VonMisesStress);
title('Stress distribution')
colormap('jet')
%% Multisphere example IV (Compression)
close all; clc; 

start = tic;
[Ri, Ro] = deal(.10, .19);  % 50cm 70 cm

% material moduli for internal and external IAB skin
[C1, C2] = deal(11e11, 22e9);
% IAB material density
[rho, nu] = deal(2e-5, .495);  % 45/386 cue from dynamics of recatangular block example
ri = .08; % meters
close all
[P, model, ro, result] = bvp_free(C1, C2, Ri, Ro, rho, nu, 'compress', ri)
fprintf('Time to run: %f seconds', toc(start))

%homedir = '/Users/olalekanogunmolu/';
full_path=fullfile(homedir, savedir);
% Plot the displacements
figure(1)
subplot(131)
pdeplot3D(model,'ColorMapData',result.Displacement.ux*10);
title('x-displacement')
% : R_i = 50cm, R_o = 70cm, r_i = R_i + 20 (cm)
colormap('jet')
% colorbar('off')

subplot(132)
pdeplot3D(model,'ColorMapData',result.Displacement.uy*10);
xlabel('Time')
title('y-displacement')
colormap('jet')
% colorbar('off')

subplot(133)
pdeplot3D(model,'ColorMapData',result.Displacement.uz*10);
title('z-displacement')
% : R_i = 50cm, R_o = 70cm, r_i = R_i + 20 (cm)
colormap('jet')

figure(2)
subplot(131)
pdegplot(model,'CellLabels','on','FaceAlpha',0.6,'FaceLabels','off')
title('IAB Skins');

subplot(132)
pdeplot3D(model)
title('Mesh Model')

subplot(133)
pdeplot3D(model,'ColorMapData',result.VonMisesStress);
title('Stress distribution')
colormap('jet')