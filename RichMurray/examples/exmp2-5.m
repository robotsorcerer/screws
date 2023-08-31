(*
 * Example 2.2:  Twist coordinates for rotation about an axis
 *
 * RMM, 4 Nov 93
 *
 *)

<<Screws.m

(* Define some trig simplification rules *)
trigexpand = { Sin[x_]^2 -> (1 - Cos[x]^2) };
trigreduce = { (1 - Cos[x_]^2) -> Sin[x]^2 };
trigSimplify[expr_] := Simplify[expr /. trigexpand, Trig->False] /. trigreduce;

(*
 * First part: compute the twist coords for absolute transformation
 *
 *)

(* Problem set up: construct rigid body transformation *)
w = {0,0,1};
Rab = SkewExp[w, th[t]];
pab = {-l2 Sin[th[t]], l1 + l2 Cos[th[t]], l0};
gab = RPToHomogeneous[Rab, pab];

(* Figure out the velocities *)
Vshat = trigSimplify[ D[gab, t]. RigidInverse[gab] ];
Vs = HomogeneousToTwist[Vshat]

Vbhat = trigSimplify[ RigidInverse[gab] . D[gab, t]];
Vb = HomogeneousToTwist[Vbhat]

