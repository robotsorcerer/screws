(*
 * Example 4.3:  Dynamics of a three-link manipulator
 *
 * RMM, 5 Nov 93
 *
 *)

<<RobotLinks.m

(* Define configuration and velocity vectors *)
q = {th1,th2,th3}
w = {dth1,dth2,dth3}

(* Define the twists which define the kinematics *)
xi1 = RevoluteTwist[{0,0,0}, {0,0,1}];
xi2 = RevoluteTwist[{0,0,l0}, {-1,0,0}];
xi3 = RevoluteTwist[{0,l1,l0}, {-1,0,0}];
xi0 = {0,0,0, 0,0,0};		(* zero twist *)

(* Compute the kinematics for the link frames *)
g10 = RPToHomogeneous[IdentityMatrix[3], {0,0,r0}];
g1 = ForwardKinematics[{xi1,th1}, g10];
J1 = Simplify[ BodyJacobian[{xi1,th1}, {xi0,th2}, {xi0,th3}, g10] ];

g20 = RPToHomogeneous[IdentityMatrix[3], {0,r1,l0}];
g2 = ForwardKinematics[{xi1,th1}, {xi2,th2}, g20];
J2 = Simplify[ BodyJacobian[{xi1,th1}, {xi2,th2}, {xi0,th3}, g20] ];

g30 = RPToHomogeneous[IdentityMatrix[3], {0,l1+r2,l0}];
g3 = ForwardKinematics[{xi1,th1}, {xi2,th2}, {xi3,th3}, g30];
J3 = Simplify[ BodyJacobian[{xi1,th1}, {xi2,th2}, {xi3,th3}, g30] ];

(* Now compute the inertia matrix for the manipulator *)
M1 = DiagonalMatrix[{m1,m1,m1,Ix1,Iy1,Iz1}];
M2 = DiagonalMatrix[{m2,m2,m2,Ix2,Iy2,Iz2}];
M3 = DiagonalMatrix[{m3,m3,m3,Ix3,Iy3,Iz3}];

Inertia = Simplify[
  Transpose[J1].M1.J1 + Transpose[J2].M2.J2 + Transpose[J3].M3.J3
]

(* Compute the Coriolis matrix using RobotLinks *)
Coriolis = InertiaToCoriolis[Inertia, q, w];

(* Extract back out the Christoffel symbols *)
gamma[i_,j_,k_] := D[Coriolis, w[[k]]][[i,j]]

(* Print out all of the Christoffel symbols *)
Table[
  Print[
    "gamma[", i, ",", j, ",", k, "]  =",
    InputForm[Simplify[gamma[i,j,k], Trig->False]]
  ],
  {i,1,3}, {j,1,3}, {k,1,3}
];

(* Compute the forces due to gravity *)
h1 = Simplify[ RigidPosition[g1][[3]] ];
h2 = Simplify[ RigidPosition[g2][[3]] ];
h3 = Simplify[ RigidPosition[g3][[3]] ];

V = Simplify[ m1 g h1 + m2 g h2 + m3 g h3]
Gravity = Map[D[V, #]&, q];
