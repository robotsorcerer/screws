(*
 * Example 2.6: Velocity of a two-link mechanism
 *
 * RMM, 4 Nov 93
 *
 *)

<<Screws.m

(* Define some trig simplification rules *)
trigexpand = { Sin[x_]^2 -> (1 - Cos[x]^2) };
trigreduce = { (1 - Cos[x_]^2) -> Sin[x]^2 };
trigSimplify[expr_] := Simplify[expr /. trigexpand, Trig->False] /. trigreduce;

xiab = {0,0,0, 0,0,1};  Vab = xiab dth1;
xibc = {l1,0,0, 0,0,1}; Vbc = xibc dth2;

gab = TwistExp[xiab, th1] . RPToHomogeneous[IdentityMatrix[3], {0,0,l0}];

Vac = Vab + RigidAdjoint[gab] . Vbc
