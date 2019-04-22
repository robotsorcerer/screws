(*
 * Example 4.4:  Dynamics of an idealized SCARA manipulator
 *
 * Richard M. Murray
 * 22 January 1992
 *
 * This file contains a portion of the file scara.m
 *)

<<RobotLinks.m

(* define configuration and velocity variabales *)
q = {th1, th2, th3, th4};
w = {dth1, dth2, dth3, dth4}

(* twist axes for SCARA robot, reference frame at base [from exmp3-1.m] *)
xi0 = {0,0,0, 0,0,0};		(* dummy twist *)
xi1 = {0,0,0, 0,0,1};		(* base *)
xi2 = {l1,0,0, 0,0,1};		(* elbow *)
xi3 = {l1+l2,0,0, 0,0,1};	(* wrist revolute *)
xi4 = {0,0,1, 0,0,0};		(* wrist prismatic *)

(* Define the link frames and compute body Jacobians *)
g10 = RPToHomogeneous[IdentityMatrix[3], {0,r1,l0}];
J1 = Simplify[ BodyJacobian[{xi1,th1}, {xi0,th2}, {xi0,th3}, {xi0,th4}, g10] ];

g20 = RPToHomogeneous[IdentityMatrix[3], {0,l1+r2,l0+h2}];
J2 = Simplify[ BodyJacobian[{xi1,th1}, {xi2,th2}, {xi0,th3}, {xi0,th4}, g20] ];

g30 = RPToHomogeneous[IdentityMatrix[3], {0,l1+l2,l0+h3}];
J3 = Simplify[ BodyJacobian[{xi1,th1}, {xi2,th2}, {xi3,th3}, {xi0,th4}, g30] ];

g40 = RPToHomogeneous[IdentityMatrix[3], {0,l1+l2,l0+h3+h4}];
J4 = Simplify[ BodyJacobian[{xi1,th1}, {xi2,th2}, {xi3,th3}, {xi4,th4}, g40] ];

(* Compute the inertia matrix *)
M1 = DiagonalMatrix[{m1,m1,m1,Ix1,Iy1,Iz1}];
M2 = DiagonalMatrix[{m2,m2,m2,Ix2,Iy2,Iz2}];
M3 = DiagonalMatrix[{m3,m3,m3,Ix3,Iy3,Iz3}];
M4 = DiagonalMatrix[{m4,m4,m4,Ix4,Iy4,Iz4}];

(* Put the matrix together and make paramter substitutions *)
Inertia = Simplify[
  Transpose[J1].M1.J1 + Transpose[J2].M2.J2 + 
  Transpose[J3].M3.J3 + Transpose[J4].M4.J4
  /. {
    Iz1 -> alpha - r1^2 m1 - l1^2 (m2 + m3 + m4) ,
    Iz2 -> beta - delta - l2^2 (m3 + m4) - m2 r2^2,
    x___ l1 l2 m3 y___ -> x (gamma - l1 l2 m4 - l1 m2 r2) y,
    Iz3 -> delta - Iz4
  }
];

(* Compute the Coriolis matrix using RobotLinks *)
Coriolis = InertiaToCoriolis[Inertia, q, w];

