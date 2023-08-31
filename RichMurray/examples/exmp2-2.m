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
Rab = SkewExp[w, alpha];
pab = {-l2 Sin[alpha], l1 + l2 Cos[alpha], 0};
gab = RPToHomogeneous[Rab, pab];

(* Figure out the vector v *)
A = (IdentityMatrix[3] - Rab) . Skew[w] + Transpose[{w}] . {w} alpha;
Ainv = trigSimplify[ Inverse[A] ];
v = trigSimplify[Ainv . pab];

(* Compute the rigid motion associated with a twist *)
xi = {l1,0,0,0,0,1};
expxi = TwistExp[xi, th]
gab0 = RPToHomogeneous[IdentityMatrix[3], {0,l1,0}];

(* figure out the twist coordinates *)
xi = Join[v, w]

(*
 * Second part: compute the twist coords for relative transformation
 *
 *)

grel = trigSimplify[gab . RigidInverse[gab /. alpha-> 0]];
xirel = RigidTwist[grel /. alpha -> 1]

(*! This generates a strange error message which I don't understand. -rmm !*)
