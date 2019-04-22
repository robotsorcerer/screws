(*
 * Example 3.17:  Self-motion manifold for a planar manipulator
 *
 * RMM, 21 Sep 92
 *
 *)

<< RobotLinks.m

(*Compute the forward kinematics map *)
gst = Simplify[
  TwistExp[RevoluteTwist[{0,0,0}, {0,0,1}], th1] .
  TwistExp[RevoluteTwist[{l1,0,0}, {0,0,1}], th2] .
  TwistExp[RevoluteTwist[{l1+l2,0,0}, {0,0,1}], th3] .
  RPToHomogeneous[IdentityMatrix[3], {l1+l2+l3,0,0}]
];

(* Now extract out the xy position of the end-effector *)
pst = RigidPosition[gst];
p = pst[[{1,2}]];  

(*
 * Compute the Jacobian of the kinematic map 
 *
 * Note: this is *not* the manipulator Jacobian since we are not
 * modelling the forward kinematics as a map g:Q -> SE(3).
 *
 *)

(* Utility function for computing the Jacobian of a mapping *)
Jac[f_, x_] := 
  If[VectorQ[ f ],
    Table[ D[ f[[i]], x[[j]] ], {i, Length[f]} , {j, Length[x]} ],
  (* else *)
    Table[ D[ f, x[[j]] ], {j, Length[x]} ]
  ];

(* Now calculate the Jacobian and its null space *)
J = Jac[p, {th1,th2,th3}];
vN = Simplify[NullSpace[J]]


(* 
 * Plot out the self-motion manifold
 *
 * )

(* Define a couple of utility functions *)
mod[x_] := N[x - Floor[(x+Pi)/(2Pi)] * 2Pi];
mod1[x_] := N[x - Floor[x/(2Pi)] * 2Pi];

(* Function to solve the inverse kinematics of a two link manipulator *)
inverse[x_, y_, flip_] :=
  Module[
    {th1, th2, alpha, r = Sqrt[x^2 + y^2]},
    alpha = ArcCos[(2 - r^2)/2];    th2 = Pi + flip alpha;
    beta = ArcCos[r/2];		    th1 = ArcTan[y/x] + flip beta;
    {mod[N[th1]], mod[N[th2]]}
  ];

(* Now swing the last joint around and solve for remaining variables *)
genpt[th3_] :=
  Join[
    inverse[1 + Sin[th3], 1 + Cos[th3], If[N[th3] > N[2Pi], -1, 1]],
    {th3}]

ParametricPlot3D[{genpt[th][[1]], genpt[th][[2]], mod1[th]}, {th, 0, 4Pi}, PlotPoints->100]
