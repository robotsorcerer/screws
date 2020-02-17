%
%PGD Code for hyperelasticity
% I. Alfaro , D. Gonzalez , E. Cueto
% Universidad de Zaragoza
% AMB - I3A Dec 2015
%
clear all ; clc ; format long g; close all ;
%
%VARIABLES
%%
%Global variables .
global coords tet tri tf
global Fprev Gprev num_iter_prev
global E nu
global s1 vx vp
E = 2.1e11; nu = 0.25; % Material ( Young Modulus and Poisson Coef ).
Modulus_init = 30e6;% Total force Modulus .
nincr = 5; % Number of load increments .
Modulus = Modulus_init / nincr ; % Force on each load increment .
TOL = 0.05; % Tolerance .
iter = zeros(1); % # of iterations needed on enrichment .m function

num_max_iter = 10; % Max . # of summands for the approach on each load incr .
%
%GEOMETRY AND BOUNDARY CONDITIONS
%
%Nodes and elements readed from external files .
coords = load('gcoordBeam.dat'); % Nodal coordinates .
tet = load (’conecBeam .dat ’); % Connectivity list .
Ind = 1: size (coords ,1); % List of nodes .
bcnode = Ind ( coords (: ,1)== min ( coords (: ,1))); % Boundary : Fix left side .
IndBcnode = sort ([3*( bcnode -1)+1 3*( bcnode -1)+2 3* bcnode ]); % D.o.f. BCs .
% Load can be applied in every node of the surface
% Make use of triangulation MatLab function to obtain boundary surface .
TR = triangulation (tet , coords );
% Connectivity of tf corresponds to the node number of the whole domain
% Connectivity of tri corresponds to the node number of the free boundary
[tf] = freeBoundary (TR ); % Dependent of 3D geometry of the boundary .
[tri , coors ] = freeBoundary (TR ); % Independent triangulation of boundary .
IndS = 1: size (coors ,1); ncoors = numel ( IndS );
%%
ALLOCATION OF MATRICES AND VECTORS
%%
Vector solutions of all previous increments ( cumulative ).
Fprev = zeros ( numel ( coords ) ,1); % Vector solution for space .
Gprev = zeros ( size (coors ,1) ,1); % Vector solution for load .
num_iter_prev = 0; % total # of summands .
% Loop on every load increment
for incr =1: nincr
fprintf (1, ’load ␣ step ␣%d\n’,incr );
%%
INITIALIZATION OF MATRICES AND VECTORS
%
F = zeros ( numel ( coords ) ,1);% For space in current increment .
R = zeros ( numel ( coords ) ,1);% For load in current increment .
G = zeros ( size (coors ,1) ,1);
S = zeros ( size (coors ,1) ,1);
num_iter = 0;
Error_iter = 1.0;
Aprt = 0;
%%
STIFNESS AND MASS MATRICES COMPUTATION
%
fem3DHyperelastic ;
coorp = 1: size (coors ,1); % We consider each position like cases of load
elemstiffHyperelastic ( coorp );
%%
SOURCE
%%
Identifying local nodes of the loaded surface on the global connectivity .
% To obtain that : coords (IndL ,:) - coors = zeros (nn2 ,1).
[trash ,trash2 ,xj] = intersect (IndS ,tri (:)); % TRI Local connectivity .
IndL = tf(xj ); % TF Global connectivity of the loaded surface .
DOFLoaded = 3* IndL ; % Consider vertical load on the top of the beam .
vx = zeros ( numel ( coords ), ncoors );
vx( DOFLoaded ,:) = -Modulus .* eye ( ncoors ); % Space terms for the source .