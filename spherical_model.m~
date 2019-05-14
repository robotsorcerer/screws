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

%% Deformation analysis
clc; clear all; close all
homedir = '/home/lex/';
otherpath = 'Documents/superchicko/sofa/IAB8/compoz';
sphere_path=fullfile(homedir, otherpath, 'sphere.STL');
% model = createpde(1);
model = createpde('structural', 'static-solid');
importGeometry(model, sphere_path);
%% Plots
close all
figure(1)
subplot(211)
pdegplot(model,'FaceLabels','on')
view(30,30);
title('Frontal view')

subplot(212)
pdegplot(model,'FaceLabels','on')
view(30,30);
title('Rear view')

% specify material moduli, poison ratio for the deformation;
structuralProperties(model,'YoungsModulus',200e9,...
            'PoissonsRatio',.45,'MassDensity',.5/386);
% define the boundary conditions
structuralBC(model,'Face',1,'Constraint','symmetric'); % symmetric isochoric

%Specify the gravity load on the sphere.
structuralBodyLoad(model,'GravitationalAcceleration',[0;0;-9.8]);

%% Multisphere; Concentric sphere with internal cavity
close all; clc
gm = multisphere([50 67])
model = createpde('structural','static-solid')
model.Geometry = gm

pdegplot(model,'CellLabels','on','FaceAlpha',0.4,'FaceLabels','on')
title('IAB with incompressibility constraints');
% For each SoRo, specify the Young's modulus, Poisson ratio and mass
% density
structuralProperties(model, 'Cell', 1, 'YoungsModulus',110E9,...
            'PoissonsRatio',.45,'MassDensity',45);
structuralProperties(model, 'Cell', 2, 'YoungsModulus',200E9,...
            'PoissonsRatio',.45,'MassDensity',45);
% define the boundary conditions
structuralBC(model,'Face',1,'Constraint','symmetric'); % symmetric isochoric
structuralBC(model,'Face',2,'Constraint','symmetric'); % symmetric isochoric

% Specify boundary loads for structural model
%Specify the gravity load on the sphere.
structuralBodyLoad(model, 'GravitationalAcceleration',[0;0;-9.8]);

% Specify surface traction for faces 1 and 2
structuralBoundaryLoad(model,'Face',[1,2],'SurfaceTraction',[0;0;100])

%specify contact-free pressure as a function handle (see thesis)
int_pressure = @(location,state)5E7.*sin(25.*state.time);

structuralBoundaryLoad(model,'Face',2,'Pressure',int_pressure)
